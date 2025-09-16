#define SDL_MAIN_HANDLED
#include <SDL2/SDL.h>
#include <stdbool.h>

#include "res/Brickwall3_Texture.h"

uint32_t *textures[1] = {
    Brickwall3_Texture};

SDL_Window *window;
SDL_Surface *surface;

uint32_t *pixels;
bool keys[512];

// uint32_t rgb(uint8_t r, uint8_t g, uint8_t b) {
// 	return (uint32_t)SDL_MapRGB(surface->format, r, g, b);
// }

unsigned short scrW = 640;
unsigned short scrH = 320;

#define FOCAL 200.0f
#define PI3 M_PI *M_PI *M_PI

void pixel(int x, int y, uint32_t color)
{
    y = -y + scrH;
    if (x < 0 || x >= scrW || y < 0 || y >= scrH)
        return;
    pixels[y * scrW + x] = color;
}

typedef struct
{
    double x, y, z, a, l;
    double sa, ca;
} player_t;

player_t p;

typedef struct
{
    double x1, y1,
        x2, y2;
    int c;
} wall_t;

typedef struct
{
    int idx, end;
    double z1, z2;
    int c;
    double x, y;
    double dist;
    int *slut;
    signed char surf;
} sector_t;

wall_t *walls;
sector_t *sectors;
int wallCount = 0;
int sectCount = 0;

int loadMap(const char *filePath)
{
    FILE *f = fopen(filePath, "r");
    char ln[256];
    while (fgets(ln, 256, f))
    {
        switch (ln[0])
        {
        case 'w':
            wallCount += 1;
            walls = realloc(walls, sizeof(wall_t) * (wallCount));
            walls[wallCount - 1] = (wall_t){0};
            sscanf(ln,
                   "w %lf %lf %lf %lf %d",
                   &walls[wallCount - 1].x1,
                   &walls[wallCount - 1].y1,
                   &walls[wallCount - 1].x2,
                   &walls[wallCount - 1].y2,
                   &walls[wallCount - 1].c);
            break;
        case 's':
            sectCount += 1;
            sectors = realloc(sectors, sizeof(sector_t) * (sectCount));
            sectors[sectCount - 1] = (sector_t){0};
            sscanf(ln,
                   "s %d %d %lf %lf %d",
                   &sectors[sectCount - 1].idx,
                   &sectors[sectCount - 1].end,
                   &sectors[sectCount - 1].z1,
                   &sectors[sectCount - 1].z2,
                   &sectors[sectCount - 1].c);
            break;
        case 'p':
            p = (player_t){0};
            sscanf(ln, "p %lf %lf %lf %lf %lf", &p.x, &p.y, &p.z, &p.a, &p.l);
            break;
        default:
            break;
        }
        printf(ln);
    }
    fclose(f);
    return 1;
}

double dist(double x, double y)
{
    return sqrt((x * x) + (y * y));
}

// void floors()
// {
//     int x, y;

// }

void drawWall(int x1, int x2, int b1, int b2, int t1, int t2, int s, int w, int i)
{
    int x, y;

    float ht = 0;
    float ht_step = (float)BRICKWALL3_TEXTURE_WIDTH / (float)(x2 - x1);

    // dist
    int dyb = b2 - b1;
    int dyt = t2 - t1;
    int dx = x2 - x1;
    if (!dx)
        dx = 1;

    int xs = x1; // x start

    if (x1 < 0)
    {
        ht -= ht_step * x1;
    }
    x1 = SDL_clamp(x1, 0, scrW - 1);
    x2 = SDL_clamp(x2, 0, scrW - 1);

    for (x = x1; x < x2; x += 1)
    {
        int y1 = dyb * (x - xs + 0.5) / dx + b1;
        int y2 = dyt * (x - xs + 0.5) / dx + t1;

        float vt = 0;
        float vt_step = (float)BRICKWALL3_TEXTURE_HEIGHT / (float)(y2 - y1);

        if (y1 < 0)
        {
            vt -= vt_step * y1;
        }
        y1 = SDL_clamp(y1, 0, scrH - 1);
        y2 = SDL_clamp(y2, 0, scrH - 1);

        if (i == 0)
        {
            if (sectors[s].surf == 1)
                sectors[s].slut[x] = y1;
            if (sectors[s].surf == 2)
                sectors[s].slut[x] = y2;

            for (y = y1; y < y2; y += 1)
            {
                pixel(x, y, Brickwall3_Texture[(int)vt * 64 + (int)ht]);
                vt += vt_step;
            }
        }
        if (i == 1)
        {
            int xo = scrW / 2;
            int yo = scrH / 2;

            int x2 = x - xo;
            int wo;

            if (sectors[s].surf == 1) {
                y2 = sectors[s].slut[x];
                wo = sectors[s].z1;
            }
            if (sectors[s].surf == 2) {
                y1 = sectors[s].slut[x];
                wo = sectors[s].z2;
            }

            float la = p.l * (M_PI * 2);
            if (la > scrH)
                la = scrH;
            float za = (p.z - wo) / yo;
            if (!za)
                za = 0.0001;

            int ystart = y1 - yo;
            int yend = y2 - yo;

            for (y = ystart; y < yend; y += 1)
            {
                    float z = y + la;
                    if (!z)
                        z = 0.0001;
                    float fx = x2 / z * za;
                    float fy = FOCAL / z * za;

                    float rx = fx * p.ca + fy * p.sa - (p.x / (120.0));
                    float ry = fy * p.ca - fx * p.sa - (p.y / (120.0));

                    if (rx < 0) rx = -rx + 1;
                    if (ry < 0)  ry = -ry + 1;

                    pixel(x2 + xo, y + yo, Brickwall3_Texture[((int)ry % 64) + ((int)rx % 64)]);
            }
        }
        ht += ht_step;
    }
}

void clipBehind(int *x1, int *y1, int *z1, int x2, int y2, int z2)
{
    double da = *y1;
    double d = da - y2;
    if (!d)
        d = 1;
    double s = da / (da - y2);

    *x1 = *x1 + s * (x2 - (*x1));
    *y1 = *y1 + s * (y2 - (*y1));
    if (!(*y1))
        *y1 = 1;
    *z1 = *z1 + s * (z2 - (*z1));
}

void render()
{
    int s = 0, w;

    double x = sectors[s].x - p.x;
    double y = sectors[s].y - p.y;

    double oldx = x;
    x = x * p.ca + y * p.sa;
    y = y * p.ca - oldx * p.sa;

    sectors[s].dist = dist(x, y);

    for (s = 1; s != sectCount; s += 1)
    {
        x = sectors[s].x - p.x;
        y = sectors[s].y - p.y;

        oldx = x;
        x = x * p.ca + y * p.sa;
        y = y * p.ca - oldx * p.sa;

        sectors[s].dist = dist(x, y);
        if (sectors[s].dist > sectors[s - 1].dist)
        {
            sector_t swp;
            swp = sectors[s];
            sectors[s] = sectors[s - 1];
            sectors[s - 1] = swp;
        }
    }
    for (s = 0; s != sectCount; s += 1)
    {
        int cycles;
        if (p.z < sectors[s].z1)
        {
            sectors[s].surf = 1;
            cycles = 2;
            for (int x = 0; x < scrW; x += 1)
                sectors[s].slut[x] = scrH;
            // memset(sectors[s].slut, scrH, sizeof(int) * scrW);
        }
        else if (p.z > sectors[s].z2)
        {
            sectors[s].surf = 2;
            cycles = 2;
            for (int x = 0; x < scrW; x += 1)
                sectors[s].slut[x] = 0;
            //            memset(sectors[s].slut, 0, sizeof(int) * scrW);
        }
        else
        {
            sectors[s].surf = 0;
            cycles = 1;
        }

        int wx[4], wy[4], wz[4];
        for (int i = 0; i != cycles; i += 1)
        {
            for (w = sectors[s].idx; w != sectors[s].end; w += 1)
            {
                // translation
                int x1 = walls[w].x1 - p.x,
                    y1 = walls[w].y1 - p.y;

                int x2 = walls[w].x2 - p.x,
                    y2 = walls[w].y2 - p.y;

                // on i = 1, we flip the faces, doing the backface
                if (i)
                {
                    int swap = x1;
                    x1 = x2;
                    x2 = swap;

                    swap = y1;
                    y1 = y2;
                    y2 = swap;
                }

                // rotate
                wx[0] = x1 * p.ca - y1 * p.sa;
                wx[1] = x2 * p.ca - y2 * p.sa;
                wx[2] = wx[0];
                wx[3] = wx[1];

                wy[0] = y1 * p.ca + x1 * p.sa;
                wy[1] = y2 * p.ca + x2 * p.sa;
                wy[2] = wy[0];
                wy[3] = wy[1];

                wz[0] = sectors[s].z1 - p.z - p.l * wy[0] / 32.0; // / 32.0;
                wz[1] = sectors[s].z1 - p.z - p.l * wy[1] / 32.0; // / 32.0;
                wz[2] = sectors[s].z2 - p.z - p.l * wy[2] / 32.0; // / 32.0;
                wz[3] = sectors[s].z2 - p.z - p.l * wy[3] / 32.0; // / 32.0;

                // clip
                if (wy[0] <= 0 && wy[1] <= 0)
                    continue;

                if (wy[0] <= 0)
                {
                    clipBehind(&wx[0], &wy[0], &wz[0], wx[1], wy[1], wz[1]);
                    clipBehind(&wx[2], &wy[2], &wz[2], wx[3], wy[3], wz[3]);
                }
                if (wy[1] <= 0)
                {
                    clipBehind(&wx[1], &wy[1], &wz[1], wx[0], wy[0], wz[0]);
                    clipBehind(&wx[3], &wy[3], &wz[3], wx[2], wy[2], wz[2]);
                }

                // screen space

                wx[0] = wx[0] * FOCAL / wy[0] + scrW / 2;
                wx[1] = wx[1] * FOCAL / wy[1] + scrW / 2;

                wy[0] = wz[0] * FOCAL / wy[0] + scrH / 2;
                wy[1] = wz[1] * FOCAL / wy[1] + scrH / 2;
                wy[2] = wz[2] * FOCAL / wy[2] + scrH / 2;
                wy[3] = wz[3] * FOCAL / wy[3] + scrH / 2;

                drawWall(wx[0], wx[1], wy[0], wy[1], wy[2], wy[3], s, w, i);
            }
            // sectors[s].surf *= -1;
        }
    }
}

int main()
{
    printf("%d", loadMap("C:\\Users\\User\\Pictures\\doom-style-game\\map.txt"));

    SDL_Init(SDL_INIT_VIDEO);

    SDL_Rect Bounds;
    SDL_GetDisplayBounds(0, &Bounds);
    scrW = Bounds.w;
    scrH = Bounds.h;

    window = SDL_CreateWindow("doom-style game!", 0, 0, scrW, scrH, SDL_WINDOW_BORDERLESS);
    surface = SDL_GetWindowSurface(window);
    pixels = (uint32_t *)surface->pixels;

    for (int s = 0; s != sectCount; s += 1)
    {
        //("%d %d %d %d", sectors[s].idx, sectors[s].end, sectors[s].z1, sectors[s].z2);
        // calculates the sector position as a vector
        double vx;
        double vy;
        for (int w = sectors[s].idx; w != sectors[s].end; w += 1)
        {
            vx += (walls[w].x2 + walls[w].x1) / 2;
            vy += (walls[w].y2 + walls[w].y1) / 2;
        }
        sectors[s].x = vx / (sectors[s].end - sectors[s].idx);
        sectors[s].y = vy / (sectors[s].end - sectors[s].idx);

        sectors[s].slut = malloc(sizeof(int) * scrW);
    }

    double deltaTime, lastFrameTime = SDL_GetTicks(), currentTime;
#define P_SPD 0.2
#define P_RSPD 0.002
    bool quit = false;
    SDL_Event e;
    while (!quit)
    {
        currentTime = SDL_GetTicks();
        deltaTime = currentTime - lastFrameTime;
        lastFrameTime = currentTime;
        // printf("%f\n", deltaTime);
        while (SDL_PollEvent(&e))
        {
            if (e.type == SDL_QUIT)
            {
                quit = true;
            }
            else if (e.type == SDL_KEYDOWN)
            {
                keys[e.key.keysym.scancode] = true;
            }
            else if (e.type == SDL_KEYUP)
            {
                keys[e.key.keysym.scancode] = false;
            }
        }

        p.sa = sin(p.a);
        p.ca = cos(p.a);
        //
        double movex = 0, movey = 0, movez = 0;
        if (keys[SDL_SCANCODE_W])
            movey += 1;
        if (keys[SDL_SCANCODE_S])
            movey -= 1;
        if (keys[SDL_SCANCODE_A])
            movex -= 1;
        if (keys[SDL_SCANCODE_D])
            movex += 1;
        if (keys[SDL_SCANCODE_SPACE])
            movez += 1;
        if (keys[SDL_SCANCODE_LSHIFT])
            movez -= 1;
        if (keys[SDL_SCANCODE_LEFT])
            p.a -= P_RSPD * deltaTime;
        if (keys[SDL_SCANCODE_RIGHT])
            p.a += P_RSPD * deltaTime;
        if (keys[SDL_SCANCODE_UP])
            p.l += P_RSPD * (M_PI*M_PI*M_PI) * deltaTime;
        if (keys[SDL_SCANCODE_DOWN])
            p.l -= P_RSPD * (M_PI*M_PI*M_PI) * deltaTime;

        int oldx = movex;
        movex = movex * p.ca + movey * p.sa;
        movey = movey * p.ca - oldx * p.sa;

        p.x += movex * P_SPD * deltaTime;
        p.y += movey * P_SPD * deltaTime;
        p.z += movez * P_SPD * deltaTime;

        SDL_LockSurface(surface);
        memset(pixels, 0, sizeof(uint32_t) * scrW * scrH);
        render();
        SDL_UnlockSurface(surface);
        SDL_UpdateWindowSurface(window);
    }
    for (int s = 0; s != sectCount; s += 1)
    {
        free(sectors[s].slut);
    }
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}