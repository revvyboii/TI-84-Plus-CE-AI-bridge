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
 *   CLEAR       clear line, cancel, or back to the menu
 *   Physics 1 solver: pick a topic, type numbers, LEFT/RIGHT picks what to find,
 *   ENTER solves. ^ or , types E for powers of ten (3E8).
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
enum { ST_COMPOSE, ST_WAIT, ST_READ, ST_MENU, ST_TOPICS, ST_FORM };
enum { M_UPPER, M_LOWER, M_NUM };

static uint8_t state = ST_MENU;
static uint8_t home_state = ST_COMPOSE; /* where CLEAR goes back to from an answer */
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

    draw_frame("ENTER send  Y= new chat  CLEAR menu");

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

/* Send `wire` to Claude; show `echo` above the answer. */
static bool send_text(const char *wire, uint16_t wlen, const char *echo, uint16_t elen) {
    if (wlen == 0) return false;
    if (!srl_open)     { notice = "Plug the USB cable into your laptop"; dirty = true; return false; }
    if (!bridge_ready) { notice = "Bridge not running on laptop";      dirty = true; return false; }

    if (elen > 600) elen = 600;
    text[0] = '>';
    text[1] = ' ';
    memcpy(text + 2, echo, elen);
    text_len = (uint16_t)(2 + elen);
    text[text_len++] = '\n';
    text[text_len++] = '\n';
    reply_offset = text_len;
    wrap_text();

    if (!send_bytes(wire, wlen) || !send_bytes("\n", 1)) {
        notice = "Send failed - check the cable";
        dirty = true;
        return false;
    }
    state = ST_WAIT;
    dirty = true;
    return true;
}

static void send_question(void) {
    if (send_text(input, input_len, input, input_len)) input_len = 0;
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
/* Physics 1 solver (menus + fill-in form)                             */
/* ================================================================== */

#define MAX_FIELDS 6
#define MAX_FINDS  6
#define VAL_LEN    11

typedef struct { const char *label; const char *unit; } field_t;
typedef struct {
    const char *name;                 /* menu name */
    const char *desc;                 /* sent to Claude */
    uint8_t nfields;
    field_t fields[MAX_FIELDS];
    uint8_t nfinds;
    const char *finds[MAX_FINDS];
} topic_t;

static const topic_t topics[] = {
    { "1D motion (const accel)", "1D kinematics with constant acceleration", 5,
      { {"v0 (initial vel)", "m/s"}, {"v (final vel)", "m/s"}, {"a (accel)", "m/s^2"},
        {"t (time)", "s"}, {"dx (displacement)", "m"} },
      5, { "final velocity", "displacement", "time", "acceleration", "initial velocity" } },

    { "Projectile motion", "projectile motion, g = 9.8 m/s^2, no air resistance", 4,
      { {"v0 (launch speed)", "m/s"}, {"angle above horiz", "deg"}, {"launch height", "m"},
        {"landing height", "m"} },
      5, { "max height", "range", "time of flight", "impact speed", "time to top" } },

    { "Incline + friction", "block on an inclined plane with friction, g = 9.8 m/s^2", 5,
      { {"m (mass)", "kg"}, {"angle of incline", "deg"}, {"mu_k (kinetic)", ""},
        {"mu_s (static)", ""}, {"F pushing up ramp", "N"} },
      5, { "acceleration", "normal force", "friction force", "does it slide?", "F to move at const v" } },

    { "Pulley / Atwood", "two masses connected by a string over a massless frictionless pulley; m1 on a table (mu_k given) or hanging if mu_k is blank, m2 hanging", 3,
      { {"m1", "kg"}, {"m2 (hanging)", "kg"}, {"mu_k (m1 on table)", ""} },
      2, { "acceleration", "tension" } },

    { "Work & energy", "work-energy / conservation of energy problem, g = 9.8 m/s^2", 6,
      { {"m (mass)", "kg"}, {"h1 (start height)", "m"}, {"v1 (start speed)", "m/s"},
        {"h2 (end height)", "m"}, {"k (spring const)", "N/m"}, {"W by friction", "J"} },
      5, { "final speed", "height reached", "spring compression", "work done", "power needed" } },

    { "Momentum / collisions", "1D collision, positive = right", 4,
      { {"m1", "kg"}, {"v1 before", "m/s"}, {"m2", "kg"}, {"v2 before", "m/s"} },
      5, { "v after (stick together)", "v1' and v2' (elastic)", "KE lost", "impulse on m1", "avg force (needs t)" } },

    { "Circular motion", "uniform circular motion, g = 9.8 m/s^2", 4,
      { {"m (mass)", "kg"}, {"r (radius)", "m"}, {"v (speed)", "m/s"}, {"T (period)", "s"} },
      5, { "centripetal accel", "net (centripetal) force", "speed", "period", "max speed flat curve (mu)" } },

    { "Torque & rotation", "rotational dynamics / torque problem", 5,
      { {"F (force)", "N"}, {"r (lever arm)", "m"}, {"angle F to r", "deg"},
        {"I (moment inertia)", "kg m^2"}, {"t (time)", "s"} },
      4, { "torque", "angular accel", "angular velocity", "rotational KE" } },
};
#define NUM_TOPICS (sizeof topics / sizeof topics[0])

static const char *main_items[] = { "Physics 1 solver", "Ask Claude anything", "Quit" };

static uint8_t menu_sel = 0;
static uint8_t topic_sel = 0;
static uint8_t form_row = 0;     /* 0..nfields-1 = values, nfields = Find row */
static uint8_t find_sel = 0;
static char vals[MAX_FIELDS][VAL_LEN + 1];
static uint8_t val_len[MAX_FIELDS];
static char wire[700];

static void draw_list(const char *heading, const char *const *items, uint8_t count,
                      uint8_t sel, const char *footer) {
    draw_frame(footer);
    gfx_SetTextFGColor(C_DIM);
    gfx_PrintStringXY(heading, TEXT_X, BODY_Y);
    for (uint8_t i = 0; i < count; i++) {
        int y = BODY_Y + 16 + i * 20;
        char num[4] = { (char)('1' + i), ')', ' ', 0 };
        if (i == sel) {
            gfx_SetColor(C_ACCENT);
            gfx_FillRectangle_NoClip(TEXT_X - 2, (uint8_t)(y - 5), SCREEN_W - 2 * TEXT_X + 4, 18);
            gfx_SetTextFGColor(C_WHITE);
        } else {
            gfx_SetTextFGColor(C_TEXT);
        }
        gfx_PrintStringXY(num, TEXT_X + 4, y);
        gfx_PrintString(items[i]);
    }
}

static void draw_menu(void) {
    draw_list("Claude for TI-84", main_items, 3, menu_sel, "arrows + ENTER, or press 1-3");
}

static void draw_topics(void) {
    static const char *names[NUM_TOPICS];
    for (uint8_t i = 0; i < NUM_TOPICS; i++) names[i] = topics[i].name;
    draw_list("Physics 1 - pick a topic", names, NUM_TOPICS, topic_sel, "ENTER pick  CLEAR back");
}

static void draw_form(void) {
    const topic_t *tp = &topics[topic_sel];
    int y = BODY_Y;

    draw_frame("type numbers  ENTER solve  CLEAR back");
    gfx_SetTextFGColor(C_ACCENT);
    gfx_PrintStringXY(tp->name, TEXT_X, y);
    gfx_SetTextFGColor(C_DIM);
    gfx_PrintStringXY("leave unknowns blank", TEXT_X, y + 10);
    y += 26;

    for (uint8_t i = 0; i <= tp->nfields; i++) {
        bool on = (i == form_row);
        if (on) {
            gfx_SetColor(C_RULE);
            gfx_FillRectangle_NoClip(TEXT_X - 2, (uint8_t)(y - 3), SCREEN_W - 2 * TEXT_X + 4, 14);
        }
        if (i < tp->nfields) {
            gfx_SetTextFGColor(C_TEXT);
            gfx_PrintStringXY(tp->fields[i].label, TEXT_X, y);
            gfx_SetTextFGColor(on ? C_ACCENT : C_TEXT);
            char shown[VAL_LEN + 2];
            memcpy(shown, vals[i], val_len[i]);
            uint8_t n = val_len[i];
            if (on) shown[n++] = '_';
            else if (n == 0) shown[n++] = '-';
            shown[n] = 0;
            gfx_PrintStringXY(shown, 168, y);
            gfx_SetTextFGColor(C_DIM);
            gfx_PrintStringXY(tp->fields[i].unit, 270, y);
        } else {
            y += 6;
            if (on) {
                gfx_SetColor(C_RULE);
                gfx_FillRectangle_NoClip(TEXT_X - 2, (uint8_t)(y - 3), SCREEN_W - 2 * TEXT_X + 4, 14);
            }
            gfx_SetTextFGColor(C_ACCENT);
            gfx_PrintStringXY("FIND: < ", TEXT_X, y);
            gfx_PrintString(tp->finds[find_sel]);
            gfx_PrintString(" >");
        }
        y += 16;
    }
    gfx_SetTextFGColor(C_DIM);
    gfx_PrintStringXY("up/down: field  left/right: what to find", TEXT_X, y + 4);
}

static void open_topic(uint8_t t) {
    topic_sel = t;
    form_row = 0;
    find_sel = 0;
    memset(val_len, 0, sizeof val_len);
    state = ST_FORM;
    dirty = true;
}

static void append(char *buf, uint16_t *len, uint16_t cap, const char *s) {
    while (*s && *len < cap - 1) buf[(*len)++] = *s++;
}

static void solve_form(void) {
    const topic_t *tp = &topics[topic_sel];
    char echo[200];
    uint16_t w = 0, e = 0;
    bool any = false;

    append(wire, &w, sizeof wire, "Physics 1 problem (Young & Freedman, SI units): ");
    append(wire, &w, sizeof wire, tp->desc);
    append(wire, &w, sizeof wire, ". Given: ");
    append(echo, &e, sizeof echo, tp->name);
    append(echo, &e, sizeof echo, ": ");

    for (uint8_t i = 0; i < tp->nfields; i++) {
        if (val_len[i] == 0) continue;
        char v[VAL_LEN + 1];
        memcpy(v, vals[i], val_len[i]);
        v[val_len[i]] = 0;
        if (any) { append(wire, &w, sizeof wire, ", "); append(echo, &e, sizeof echo, ", "); }
        append(wire, &w, sizeof wire, tp->fields[i].label);
        append(wire, &w, sizeof wire, " = ");
        append(wire, &w, sizeof wire, v);
        append(wire, &w, sizeof wire, " ");
        append(wire, &w, sizeof wire, tp->fields[i].unit);
        /* short echo: label up to the first space, e.g. "v0=20" */
        for (const char *p = tp->fields[i].label; *p && *p != ' ' && e < sizeof echo - 1; p++)
            echo[e++] = *p;
        append(echo, &e, sizeof echo, "=");
        append(echo, &e, sizeof echo, v);
        any = true;
    }
    if (!any) { notice = "Enter at least one number"; dirty = true; return; }

    append(wire, &w, sizeof wire, ". Find: ");
    append(wire, &w, sizeof wire, tp->finds[find_sel]);
    append(wire, &w, sizeof wire, ". Answer format, one short line each: EQN: the equation(s) used. "
           "SUB: numbers plugged in. ANS: result with units and 3 sig figs. "
           "If something needed is missing, say what. Max 8 lines.");
    append(echo, &e, sizeof echo, " | find ");
    append(echo, &e, sizeof echo, tp->finds[find_sel]);

    home_state = ST_FORM;
    if (srl_open && bridge_ready) {
        static const char fresh[2] = { NEW_CHAT, '\n' };  /* each problem starts a fresh chat */
        send_bytes(fresh, 2);
    }
    send_text(wire, w, echo, e);
}

static void form_keys(void) {
    const topic_t *tp = &topics[topic_sel];

    if (pressed_repeat(kb_KeyDown)) {
        if (form_row < tp->nfields) form_row++;
        dirty = true;
    } else if (pressed_repeat(kb_KeyUp)) {
        if (form_row > 0) form_row--;
        dirty = true;
    } else if (PRESSED(kb_KeyRight)) {
        find_sel = (uint8_t)((find_sel + 1) % tp->nfinds);
        dirty = true;
    } else if (PRESSED(kb_KeyLeft)) {
        find_sel = (uint8_t)((find_sel + tp->nfinds - 1) % tp->nfinds);
        dirty = true;
    } else if (PRESSED(kb_KeyEnter)) {
        solve_form();
    } else if (PRESSED(kb_KeyClear)) {
        if (form_row < tp->nfields && val_len[form_row] > 0) val_len[form_row] = 0;
        else state = ST_TOPICS;
        dirty = true;
    } else if (form_row < tp->nfields) {
        uint8_t r = form_row;
        char c = 0;
        if (pressed_repeat(kb_KeyDel)) {
            if (val_len[r] > 0) { val_len[r]--; dirty = true; }
            return;
        }
        if      (PRESSED(kb_Key0)) c = '0';
        else if (PRESSED(kb_Key1)) c = '1';
        else if (PRESSED(kb_Key2)) c = '2';
        else if (PRESSED(kb_Key3)) c = '3';
        else if (PRESSED(kb_Key4)) c = '4';
        else if (PRESSED(kb_Key5)) c = '5';
        else if (PRESSED(kb_Key6)) c = '6';
        else if (PRESSED(kb_Key7)) c = '7';
        else if (PRESSED(kb_Key8)) c = '8';
        else if (PRESSED(kb_Key9)) c = '9';
        else if (PRESSED(kb_KeyDecPnt)) c = '.';
        else if (PRESSED(kb_KeyChs) || PRESSED(kb_KeySub)) c = '-';
        else if (PRESSED(kb_KeyPower) || PRESSED(kb_KeyComma)) c = 'E'; /* 3E8 = 3x10^8 */
        if (c && val_len[r] < VAL_LEN) {
            vals[r][val_len[r]++] = c;
            dirty = true;
        }
    }
}

/* returns true if the user chose Quit */
static bool menu_keys(void) {
    int8_t pick = -1;
    if (PRESSED(kb_KeyDown)) { menu_sel = (uint8_t)((menu_sel + 1) % 3); dirty = true; }
    if (PRESSED(kb_KeyUp))   { menu_sel = (uint8_t)((menu_sel + 2) % 3); dirty = true; }
    if (PRESSED(kb_KeyEnter)) pick = (int8_t)menu_sel;
    if (PRESSED(kb_Key1)) pick = 0;
    if (PRESSED(kb_Key2)) pick = 1;
    if (PRESSED(kb_Key3) || PRESSED(kb_KeyClear)) pick = 2;
    if (pick == 0) { state = ST_TOPICS; dirty = true; }
    if (pick == 1) { home_state = ST_COMPOSE; state = ST_COMPOSE; dirty = true; }
    return pick == 2;
}

static void topics_keys(void) {
    static const kb_lkey_t num_keys[] = { kb_Key1, kb_Key2, kb_Key3, kb_Key4,
                                          kb_Key5, kb_Key6, kb_Key7, kb_Key8 };
    if (pressed_repeat(kb_KeyDown)) { topic_sel = (uint8_t)((topic_sel + 1) % NUM_TOPICS); dirty = true; }
    else if (pressed_repeat(kb_KeyUp)) { topic_sel = (uint8_t)((topic_sel + NUM_TOPICS - 1) % NUM_TOPICS); dirty = true; }
    else if (PRESSED(kb_KeyEnter)) open_topic(topic_sel);
    else if (PRESSED(kb_KeyClear)) { state = ST_MENU; dirty = true; }
    else {
        for (uint8_t i = 0; i < NUM_TOPICS && i < 8; i++)
            if (PRESSED(num_keys[i])) { open_topic(i); break; }
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

        if (PRESSED(kb_KeyYequ) && (state == ST_COMPOSE || state == ST_READ || state == ST_WAIT))
            start_new_chat();

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
                else { state = ST_MENU; dirty = true; }
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
                state = home_state;
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
            } else if (PRESSED(kb_KeyEnter)) {
                state = ST_COMPOSE;            /* type a follow-up */
                dirty = true;
            } else if (PRESSED(kb_KeyClear)) {
                state = home_state;            /* back to the form or typing */
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

        case ST_MENU:
            if (menu_keys()) quit = true;
            break;

        case ST_TOPICS:
            topics_keys();
            break;

        case ST_FORM:
            form_keys();
            break;
        }

        if (dirty) {
            dirty = false;
            if (state == ST_COMPOSE)      draw_compose();
            else if (state == ST_WAIT)    draw_wait(dots);
            else if (state == ST_MENU)    draw_menu();
            else if (state == ST_TOPICS)  draw_topics();
            else if (state == ST_FORM)    draw_form();
            else                          draw_read();
            gfx_SwapDraw();
        }
    }

    usb_Cleanup();
    gfx_End();
    return 0;
}
