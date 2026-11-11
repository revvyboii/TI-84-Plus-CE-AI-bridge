/*
 * CLAUDE - chat with Claude from a TI-84 Plus CE
 *
 * The calculator acts as a USB serial device (srldrvce). A bridge program on
 * your laptop reads each question, asks Claude, and sends the answer back.
 *
 * Protocol (keep in sync with bridge.js):
 *   calc -> laptop : question text, then '\n'
 *                    0x02 '\n' = start a new chat
 *   laptop -> calc : 0x06 = "bridge is running" (sent every few seconds)
 *                    answer text (printable ASCII + '\n'), then '\0'
 *
 * Keys:
 *   ALPHA       cycle input mode  ABC -> abc -> 123
 *   ENTER       send question / go back to typing
 *   DEL         backspace
 *   CLEAR       clear line, cancel, or quit
 *   Y=          new chat (Claude forgets the conversation)
 *   UP/DOWN     scroll the answer (LEFT/RIGHT = page)
 *   In ABC/abc mode: 0 = space, . = period, (-) = ?, 3 = comma, + = '
 */

#include <srldrvce.h>
#include <usbdrvce.h>
#include <graphx.h>
#include <keypadc.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

/* ---------- layout ---------- */
#define SCREEN_W    320
#define COLS        38      /* characters per line (8px monospace) */
#define TEXT_X      6
#define TITLE_H     16
#define BODY_Y      22
#define LINE_H      10
#define BODY_LINES  20
#define FOOTER_Y    228

/* ---------- limits ---------- */
#define MAX_INPUT   240
#define MAX_TEXT    4096
#define MAX_LINES   480

/* ---------- protocol bytes ---------- */
#define MSG_END      '\0'
#define BRIDGE_HELLO 0x06
#define NEW_CHAT     0x02

/* ---------- colors (palette slots we set ourselves) ---------- */
enum { C_BG, C_TEXT, C_ACCENT, C_DIM, C_WHITE, C_RULE };
#define C_NONE 254 /* unused slot = transparent text background */

static const uint16_t palette[] = {
    gfx_RGBTo1555(250, 248, 243), /* C_BG     */
    gfx_RGBTo1555( 31,  31,  31), /* C_TEXT   */
    gfx_RGBTo1555(201, 100,  66), /* C_ACCENT */
    gfx_RGBTo1555(122, 118, 109), /* C_DIM    */
    gfx_RGBTo1555(255, 255, 255), /* C_WHITE  */
    gfx_RGBTo1555(229, 224, 214), /* C_RULE   */
};

/* ---------- serial state ---------- */
static srl_device_t srl;
static uint8_t srl_buf[1024];
static bool srl_open = false;
static usb_device_t pending_dev = NULL;
static bool bridge_ready = false;
static uint32_t last_hello = 0;

/* ---------- app state ---------- */
enum { ST_COMPOSE, ST_WAIT, ST_READ };
enum { M_UPPER, M_LOWER, M_NUM };

static uint8_t state = ST_COMPOSE;
static uint8_t mode = M_UPPER;
static bool dirty = true;
static bool discard_reply = false; /* drop the answer to a cancelled question */
static bool have_answer = false;
static const char *notice = NULL;  /* one-line message shown in the footer */

static char input[MAX_INPUT + 1];
static uint8_t input_len = 0;

static char text[MAX_TEXT];        /* "> question\n\nanswer" */
static uint16_t text_len = 0;
static uint16_t reply_offset = 0;  /* where the answer starts in text[] */

static uint16_t line_start[MAX_LINES];
static uint8_t line_len[MAX_LINES];
static uint16_t line_count = 0;
static uint16_t scroll = 0;

/* ---------- time (usbdrvce owns timer 2; read it through this) ---------- */
#define MS(x) ((uint32_t)(x) * 48000UL)
static uint32_t now(void) { return usb_GetCycleCounter(); }

/* ================================================================== */
/* USB                                                                 */
/* ================================================================== */

static void close_serial(void) {
    if (srl_open) {
        srl_Close(&srl);
        srl_open = false;
    }
    bridge_ready = false;
    dirty = true;
}

static usb_error_t handle_usb_event(usb_event_t event, void *event_data,
                                    usb_callback_data_t *callback_data) {
    usb_error_t err;
    if ((err = srl_UsbEventCallback(event, event_data, callback_data)) != USB_SUCCESS)
        return err;

    /* The laptop finished setting us up as a serial device. srl_Open must
       not be called from inside this handler, so remember the device and
       open it from the main loop. */
    if (event == USB_HOST_CONFIGURE_EVENT && !srl_open) {
        usb_device_t dev = usb_FindDevice(NULL, NULL, USB_SKIP_HUBS);
        if (dev != NULL) pending_dev = dev;
    }

    if (event == USB_DEVICE_DISCONNECTED_EVENT) {
        if (srl_open && (usb_device_t)event_data == srl.dev) close_serial();
        pending_dev = NULL;
    }
    return USB_SUCCESS;
}

static void open_pending_serial(void) {
    if (pending_dev == NULL || srl_open) return;
    srl_error_t e = srl_Open(&srl, pending_dev, srl_buf, sizeof srl_buf,
                             SRL_INTERFACE_ANY, 115200);
    pending_dev = NULL;
    srl_open = (e == SRL_SUCCESS);
    if (!srl_open) notice = "USB serial failed - replug cable";
    dirty = true;
}

/* Blocking write with a timeout. Returns false if the link went away. */
static bool send_bytes(const char *data, size_t len) {
    uint32_t start = now();
    while (len > 0) {
        if (!srl_open) return false;
        int n = srl_Write(&srl, data, len);
        if (n < 0) { close_serial(); return false; }
        data += n;
        len -= (size_t)n;
        usb_HandleEvents();
        if (now() - start > MS(3000)) return false;
    }
    return true;
}

/* ================================================================== */
/* Text layout                                                         */
/* ================================================================== */

static void wrap_text(void) {
    uint16_t i = 0;
    line_count = 0;
    while (i < text_len && line_count < MAX_LINES) {
        uint16_t start = i;
        uint16_t last_space = 0xFFFF;
        uint8_t len = 0;

        while (i < text_len && text[i] != '\n' && len < COLS) {
            if (text[i] == ' ') last_space = i;
            i++;
            len++;
        }

        if (i < text_len) {
            if (text[i] == '\n') {
                i++;                                   /* hard line break */
            } else if (len == COLS) {
                if (text[i] == ' ') {
                    i++;                               /* break right at a space */
                } else if (last_space != 0xFFFF && last_space > start) {
                    len = (uint8_t)(last_space - start);  /* word wrap */
                    i = last_space + 1;
                }                                      /* else: split a very long word */
            }
        }
        line_start[line_count] = start;
        line_len[line_count] = len;
        line_count++;
    }
}

static uint16_t max_scroll(void) {
    return line_count > BODY_LINES ? (uint16_t)(line_count - BODY_LINES) : 0;
}

/* ================================================================== */
/* Drawing                                                             */
/* ================================================================== */

static void print_n(const char *s, uint8_t n, int x, int y) {
    char tmp[COLS + 1];
    if (n > COLS) n = COLS;
    memcpy(tmp, s, n);
    tmp[n] = '\0';
    gfx_PrintStringXY(tmp, x, y);
}

static void draw_frame(const char *footer) {
    const char *status;

    gfx_FillScreen(C_BG);

    /* title bar */
    gfx_SetColor(C_ACCENT);
    gfx_FillRectangle_NoClip(0, 0, SCREEN_W, TITLE_H);
    gfx_SetTextFGColor(C_WHITE);
    gfx_PrintStringXY("Claude", TEXT_X, 4);

    if (!srl_open)          status = "no USB link";
    else if (!bridge_ready) status = "start bridge";
    else                    status = "connected";
    gfx_PrintStringXY(status, SCREEN_W - TEXT_X - (int)gfx_GetStringWidth(status), 4);

    /* footer */
    gfx_SetColor(C_RULE);
    gfx_HorizLine_NoClip(0, FOOTER_Y - 5, SCREEN_W);
    gfx_SetTextFGColor(notice ? C_ACCENT : C_DIM);
    gfx_PrintStringXY(notice ? notice : footer, TEXT_X, FOOTER_Y);
}

static void draw_compose(void) {
    static const char *mode_names[] = { "[ABC]", "[abc]", "[123]" };
    char line[COLS + 2];
    uint8_t shown = 0;
    int y = BODY_Y;

    draw_frame("ENTER send  Y= new chat  CLEAR quit");

    gfx_SetTextFGColor(C_DIM);
    gfx_PrintStringXY("Ask Claude:", TEXT_X, y);
    gfx_PrintStringXY(mode_names[mode],
                      SCREEN_W - TEXT_X - (int)gfx_GetStringWidth(mode_names[mode]), y);
    y += LINE_H + 4;

    /* input, hard-wrapped, with a cursor after the last character */
    gfx_SetTextFGColor(C_TEXT);
    for (shown = 0; shown + COLS <= input_len; shown += COLS) {
        print_n(input + shown, COLS, TEXT_X, y);
        y += LINE_H;
    }
    {
        uint8_t n = (uint8_t)(input_len - shown);
        memcpy(line, input + shown, n);
        line[n++] = '_';
        line[n] = '\0';
        gfx_PrintStringXY(line, TEXT_X, y);
        y += LINE_H;
    }

    y += LINE_H;
    gfx_SetTextFGColor(C_DIM);
    gfx_PrintStringXY("ABC/abc:  0=space  .=.  (-)=?  3=,", TEXT_X, y);
    y += LINE_H;
    if (have_answer)
        gfx_PrintStringXY("UP = see last answer", TEXT_X, y);
}

static void draw_wait(uint8_t dots) {
    static const char *spin[] = { "Thinking", "Thinking.", "Thinking..", "Thinking..." };
    draw_frame("CLEAR cancel");
    gfx_SetTextFGColor(C_ACCENT);
    for (uint16_t i = 0; i < line_count && i < BODY_LINES - 2; i++)
        print_n(text + line_start[i], line_len[i], TEXT_X, BODY_Y + (int)i * LINE_H);
    gfx_SetTextFGColor(C_DIM);
    gfx_PrintStringXY(spin[dots & 3], TEXT_X,
                      BODY_Y + (int)(line_count < BODY_LINES - 2 ? line_count + 1 : BODY_LINES - 1) * LINE_H);
}

static void draw_read(void) {
    draw_frame("arrows scroll  ENTER reply  Y= new");
    for (uint16_t row = 0; row < BODY_LINES; row++) {
        uint16_t i = scroll + row;
        if (i >= line_count) break;
        gfx_SetTextFGColor(line_start[i] < reply_offset ? C_ACCENT : C_TEXT);
        print_n(text + line_start[i], line_len[i], TEXT_X, BODY_Y + (int)row * LINE_H);
    }

    /* scrollbar */
    if (line_count > BODY_LINES) {
        int track = BODY_LINES * LINE_H;
        int bar = track * BODY_LINES / (int)line_count;
        int top = (int)((uint32_t)track * scroll / line_count);
        if (bar < 6) bar = 6;
        gfx_SetColor(C_RULE);
        gfx_FillRectangle_NoClip(SCREEN_W - 4, BODY_Y, 2, (uint8_t)track);
        gfx_SetColor(C_DIM);
        gfx_FillRectangle_NoClip(SCREEN_W - 4, (uint8_t)(BODY_Y + top), 2, (uint8_t)bar);
    }
}

/* ================================================================== */
/* Keyboard (keypadc, with edge detection and key repeat)              */
/* ================================================================== */

static uint8_t prev_keys[8];
static uint8_t new_keys[8];
static kb_lkey_t held_key = 0;
static uint32_t held_since = 0;
static uint32_t last_repeat = 0;

#define PRESSED(k) ((new_keys[(k) >> 8] & (uint8_t)(k)) != 0)
#define DOWN(k)    ((prev_keys[(k) >> 8] & (uint8_t)(k)) != 0)

static void scan_keys(void) {
    kb_Scan();
    for (uint8_t g = 1; g < 8; g++) {
        uint8_t cur = kb_Data[g];
        new_keys[g] = cur & (uint8_t)~prev_keys[g];
        prev_keys[g] = cur;
    }
}

/* true on press, then repeatedly while held */
static bool pressed_repeat(kb_lkey_t k) {
    uint32_t t = now();
    if (PRESSED(k)) {
        held_key = k;
        held_since = last_repeat = t;
        return true;
    }
    if (held_key == k) {
        if (!DOWN(k)) { held_key = 0; return false; }
        if (t - held_since > MS(400) && t - last_repeat > MS(60)) {
            last_repeat = t;
            return true;
        }
    }
    return false;
}

typedef struct { kb_lkey_t key; char alpha; char num; } keymap_t;

static const keymap_t keymap[] = {
    { kb_KeyMath,   'A',  0  }, { kb_KeyApps,  'B',  0  }, { kb_KeyPrgm,   'C',  0  },
    { kb_KeyRecip,  'D',  0  }, { kb_KeySin,   'E',  0  }, { kb_KeyCos,    'F',  0  },
    { kb_KeyTan,    'G',  0  }, { kb_KeyPower, 'H', '^' }, { kb_KeySquare, 'I',  0  },
    { kb_KeyComma,  'J', ',' }, { kb_KeyLParen,'K', '(' }, { kb_KeyRParen, 'L', ')' },
    { kb_KeyDiv,    'M', '/' }, { kb_KeyLog,   'N',  0  }, { kb_Key7,      'O', '7' },
    { kb_Key8,      'P', '8' }, { kb_Key9,     'Q', '9' }, { kb_KeyMul,    'R', '*' },
    { kb_KeyLn,     'S',  0  }, { kb_Key4,     'T', '4' }, { kb_Key5,      'U', '5' },
    { kb_Key6,      'V', '6' }, { kb_KeySub,   'W', '-' }, { kb_KeySto,    'X', '=' },
    { kb_Key1,      'Y', '1' }, { kb_Key2,     'Z', '2' }, { kb_Key3,      ',', '3' },
    { kb_Key0,      ' ', '0' }, { kb_KeyDecPnt,'.', '.' }, { kb_KeyChs,    '?', '-' },
    { kb_KeyAdd,    '\'', '+' },
};

/* Returns the character typed this frame, or 0. */
static char typed_char(void) {
    for (uint8_t i = 0; i < sizeof keymap / sizeof keymap[0]; i++) {
        if (!PRESSED(keymap[i].key)) continue;
        if (mode == M_NUM) return keymap[i].num;
        char c = keymap[i].alpha;
        if (mode == M_LOWER && c >= 'A' && c <= 'Z') c = (char)(c + 32);
        return c;
    }
    return 0;
}

/* ================================================================== */
/* Actions                                                             */
/* ================================================================== */

static void start_new_chat(void) {
    static const char msg[2] = { NEW_CHAT, '\n' };
    if (state == ST_WAIT) discard_reply = true; /* ignore the in-flight answer */
    if (srl_open) send_bytes(msg, 2);
    text_len = 0;
    line_count = 0;
    have_answer = false;
    state = ST_COMPOSE;
    notice = "New chat started";
    dirty = true;
}

static void send_question(void) {
    if (input_len == 0) return;
    if (!srl_open)     { notice = "Plug the USB cable into your laptop"; dirty = true; return; }
    if (!bridge_ready) { notice = "Bridge not running on laptop";      dirty = true; return; }

    /* echo the question at the top of the answer view */
    text[0] = '>';
    text[1] = ' ';
    memcpy(text + 2, input, input_len);
    text_len = (uint16_t)(2 + input_len);
    text[text_len++] = '\n';
    text[text_len++] = '\n';
    reply_offset = text_len;
    wrap_text();

    if (!send_bytes(input, input_len) || !send_bytes("\n", 1)) {
        notice = "Send failed - check the cable";
        dirty = true;
        return;
    }
    input_len = 0;
    state = ST_WAIT;
    dirty = true;
}

static void finish_answer(void) {
    while (text_len > reply_offset && (text[text_len - 1] == '\n' || text[text_len - 1] == ' '))
        text_len--;
    wrap_text();
    scroll = 0;
    have_answer = true;
    state = ST_READ;
    dirty = true;
}

static void poll_serial(void) {
    char buf[64];
    int n;

    if (!srl_open) return;
    n = srl_Read(&srl, buf, sizeof buf);
    if (n < 0) { close_serial(); return; }

    for (int i = 0; i < n; i++) {
        char c = buf[i];
        if (c == BRIDGE_HELLO) {
            if (!bridge_ready) dirty = true;
            bridge_ready = true;
            last_hello = now();
            continue;
        }
        if (discard_reply) {
            if (c == MSG_END) discard_reply = false;
            continue;
        }
        if (state != ST_WAIT) continue;
        if (c == MSG_END) {
            finish_answer();
        } else if (text_len < MAX_TEXT - 1) {
            text[text_len++] = c;
        }
    }

    /* bridge says hello every 2 s; if it goes quiet for 8 s, assume it stopped */
    if (bridge_ready && state != ST_WAIT && now() - last_hello > MS(8000)) {
        bridge_ready = false;
        dirty = true;
    }
}

/* ================================================================== */
/* Main                                                                */
/* ================================================================== */

int main(void) {
    uint32_t last_anim = 0;
    uint8_t dots = 0;
    bool quit = false;

    gfx_Begin();
    gfx_SetPalette(palette, sizeof palette, 0);
    gfx_SetDrawBuffer();
    gfx_SetMonospaceFont(8);
    gfx_SetTextBGColor(C_NONE);
    gfx_SetTextTransparentColor(C_NONE);

    if (usb_Init(handle_usb_event, NULL, srl_GetCDCStandardDescriptors(),
                 USB_DEFAULT_INIT_FLAGS) != USB_SUCCESS) {
        usb_Cleanup();
        gfx_FillScreen(C_BG);
        gfx_SetTextFGColor(C_TEXT);
        gfx_PrintStringXY("USB init failed. Press CLEAR.", TEXT_X, BODY_Y);
        gfx_SwapDraw();
        do kb_Scan(); while (!kb_IsDown(kb_KeyClear));
        gfx_End();
        return 1;
    }

    while (!quit) {
        usb_HandleEvents();
        open_pending_serial();
        poll_serial();
        scan_keys();

        /* any key press clears the footer notice */
        for (uint8_t g = 1; g < 8; g++)
            if (new_keys[g] && notice) { notice = NULL; dirty = true; break; }

        if (PRESSED(kb_KeyYequ)) start_new_chat();

        switch (state) {
        case ST_COMPOSE: {
            char c;
            if (PRESSED(kb_KeyAlpha)) {
                mode = (uint8_t)((mode + 1) % 3);
                dirty = true;
            } else if (PRESSED(kb_KeyEnter)) {
                send_question();
            } else if (pressed_repeat(kb_KeyDel)) {
                if (input_len > 0) { input_len--; dirty = true; }
            } else if (PRESSED(kb_KeyClear)) {
                if (input_len > 0) { input_len = 0; dirty = true; }
                else quit = true;
            } else if (PRESSED(kb_KeyUp) && have_answer) {
                state = ST_READ;
                dirty = true;
            } else if ((c = typed_char()) != 0 && input_len < MAX_INPUT) {
                input[input_len++] = c;
                dirty = true;
            }
            break;
        }

        case ST_WAIT:
            if (PRESSED(kb_KeyClear)) {
                discard_reply = true;
                text_len = 0;
                line_count = 0;
                have_answer = false;
                state = ST_COMPOSE;
                notice = "Cancelled";
                dirty = true;
            } else if (now() - last_anim > MS(350)) {
                last_anim = now();
                dots++;
                dirty = true;
            }
            break;

        case ST_READ: {
            char c;
            uint16_t ms = max_scroll();
            if (pressed_repeat(kb_KeyDown) && scroll < ms)  { scroll++; dirty = true; }
            else if (pressed_repeat(kb_KeyUp) && scroll > 0) { scroll--; dirty = true; }
            else if (PRESSED(kb_KeyRight)) {
                scroll = (uint16_t)(scroll + BODY_LINES - 1);
                if (scroll > ms) scroll = ms;
                dirty = true;
            } else if (PRESSED(kb_KeyLeft)) {
                scroll = scroll > BODY_LINES - 1 ? (uint16_t)(scroll - (BODY_LINES - 1)) : 0;
                dirty = true;
            } else if (PRESSED(kb_KeyEnter) || PRESSED(kb_KeyClear)) {
                state = ST_COMPOSE;
                dirty = true;
            } else if (PRESSED(kb_KeyAlpha)) {
                mode = (uint8_t)((mode + 1) % 3);
                state = ST_COMPOSE;
                dirty = true;
            } else if ((c = typed_char()) != 0) {
                /* start typing a follow-up right away */
                state = ST_COMPOSE;
                if (input_len < MAX_INPUT) input[input_len++] = c;
                dirty = true;
            }
            break;
        }
        }

        if (dirty) {
            dirty = false;
            if (state == ST_COMPOSE)   draw_compose();
            else if (state == ST_WAIT) draw_wait(dots);
            else                       draw_read();
            gfx_SwapDraw();
        }
    }

    usb_Cleanup();
    gfx_End();
    return 0;
}
