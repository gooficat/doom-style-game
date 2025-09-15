#define SDL_MAIN_HANDLED
#include <SDL2/SDL.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <stdio.h>


SDL_Window* window;
SDL_Surface* surface;

uint32_t* pixels;
bool keys[512];

// uint32_t rgb(uint8_t r, uint8_t g, uint8_t b) {
// 	return (uint32_t)SDL_MapRGB(surface->format, r, g, b);
// }

unsigned short scrW = 640;
unsigned short scrH = 320;


void pixel(int x, int y, uint32_t color) {
	if (x < 0 || x >= scrW || (-y+scrH) < 0 || (-y+scrH) >= scrH) return;
	pixels[(-y+scrH) * scrW + x] = color;
}

typedef struct {
    double x, y, z, a, l;
    double sa, ca;
} player_t;

player_t p;

typedef struct {
    double  x1, y1,
            x2, y2;
    uint32_t c;
} wall_t;



typedef struct {
    int idx, end;
    double z1, z2;
    uint32_t c;
    double x, y;
    double dist;
    int* slut;
    signed char surf;
} sector_t;

#define PURPLE 0xFF00FFFF
#define BLUE 0xFF00FF00
#define RED 0xFF0000FF
#define ORANGE 0xFFFF00FF

// wall_t walls[] = {
//     {0,  0, 32, 0,  PURPLE},
//     {32, 0, 32, 32, BLUE},
//     {32, 32,0,  32, PURPLE},
//     {0,  32,0,  0,  BLUE},

//     {64, 0, 96, 0,  RED},
//     {96, 0, 96, 64, ORANGE},
//     {96, 64,64, 64, RED},
//     {64, 64,64, 0,  RED},    
// };

// sector_t sectors[] = {
//     {0, 4, 0, 40, RED, BLUE},
//     {4, 8, 0, 60, ORANGE, RED},
// };

wall_t* walls;
sector_t* sectors;
int wallCount = 0;
int sectCount = 0;

int loadMap(const char* filePath) {
    FILE* f = fopen(filePath, "r");
    char ln[256];
    while (fgets(ln, 256, f)) {
        switch(ln[0]) {
            case 'w':
                wallCount++;
                walls = realloc(walls, sizeof(wall_t) * (wallCount));
                walls[wallCount-1] = (wall_t) { 0 };
                sscanf(ln,
                    "w %lf %lf %lf %lf %x",
                    &walls[wallCount-1].x1,
                    &walls[wallCount-1].y1,
                    &walls[wallCount-1].x2, 
                    &walls[wallCount-1].y2,
                    &walls[wallCount-1].c
                );
                walls[wallCount-1].c = 0xFF00FFFF;
                break;
            case 's':
                sectCount++;
                sectors = realloc(sectors, sizeof(sector_t) * (sectCount));
                sectors[sectCount-1] = (sector_t) { 0 };
                sscanf(ln,
                    "s %d %d %lf %lf %x",
                    &sectors[sectCount-1].idx,
                    &sectors[sectCount-1].end,
                    &sectors[sectCount-1].z1,
                    &sectors[sectCount-1].z2,
                    &sectors[sectCount-1].c
                );
                break;
            case 'p':
                p = (player_t) { 0 };
                sscanf(ln, "p %lf %lf %lf %lf %lf", &p.x, &p.y, &p.z, &p.a, &p.l);
                break;
            default:
                break;
        }
        // printf(ln);
    }
    fclose(f);
    return 1;
}

double dist(double x, double y) {
    // printf("%f\n", x);
    // printf("%f\n", y);
    double d = (x * x) + (y * y);
    if (!d) d = 1;
    return sqrt(d);
}

void drawWall(int x1, int x2, int b1, int b2, int t1, int t2, int s, int w, int i) {
    int x, y;

    //dist
    int dyb = b2 - b1;
    int dyt = t2 - t1;
    int dx = x2 - x1; if (!dx) dx=1;

    int xs = x1; //x start

    if (x1 < 1) x1 = 1;
    if (x1 > scrW-2) x1 = scrW - 2;
    if (x2 < 1) x2 = 1;
    if (x2 > scrW-2) x2 = scrW - 2;
    

    for (x = x1; x < x2; x++) {
        int y1 = dyb * (x - xs + 0.5) / dx + b1;
        int y2 = dyt * (x - xs + 0.5) / dx + t1;

        if (y1 < 1) y1 = 1;
        if (y1 > scrH-2) y1 = scrH - 2;
        if (y2 < 1) y2 = 1;
        if (y2 > scrH-2) y2 = scrH - 2;


        if (!i) {
            if (sectors[s].surf == 1)
                sectors[s].slut[x] = y2;
            else
                sectors[s].slut[x] = y1;
            for (y = y1; y < y2; y++)
                pixel(x, y, sectors[s].c);
        }
        else {
            if (sectors[s].surf == 1)
                y1 = sectors[s].slut[x];
            else
                y2 = sectors[s].slut[x];
            for (y = y1; y < y2; y++)
                pixel(x, y, walls[w].c);
        }
        // switch (sectors[s].surf) {
        //     case 2:
        //         sectors[s].slut[x] = y1;
        //         break;
        //     case 1:
        //         sectors[s].slut[x] = y2;
        //         break;
        //     case -2:
        //         for (y = sectors[s].slut[x]; y < y2; y++)
        //             pixel(x, y, sectors[s].cc);
        //         break;
        //     case -1: //-2
        //         for (y = y2; y < sectors[s].slut[x]; y++)
        //             pixel(x, y, sectors[s].fc);
        //         break;
        // }
//        if (sectors[s].surf < 0)

        // DrawPixel(x, y1, GREEN);
        // DrawPixel(x, y2, GREEN);
    }
}

void clipBehind(int* x1, int* y1, int* z1, int x2, int y2, int z2) {
    double da = *y1;
    double d = da - y2; if (!d) d = 1;
    double s = da / d;

    *x1 = *x1 + s * (x2 - (*x1));
    *y1 = *y1 + s * (y2 - (*y1)); if (!(*y1)) *y1 = 1;
    *z1 = *z1 + s * (z2 - (*z1));
}
#define FOCAL 250
void render() {
    int s = 0, w;

    
       double x = sectors[s].x - p.x;
       double y = sectors[s].y - p.y;
    
       double oldx = x;
        x = x * p.ca + y * p.sa;
        y = y * p.ca - oldx * p.sa;

        sectors[s].dist = dist(x, y);

    for (s = 1; s < sectCount; s++) {
        x = sectors[s].x - p.x;
        y = sectors[s].y - p.y;
    
        oldx = x;
        x = x * p.ca + y * p.sa;
        y = y * p.ca - oldx * p.sa;

        sectors[s].dist = dist(x, y);
        if (sectors[s].dist > sectors[s-1].dist) {
            sector_t swp;
            swp = sectors[s];
            sectors[s] = sectors[s-1];
            sectors[s-1] = swp;
        }
    }
    for (s = 0; s < sectCount; s++) {
//        static enum surf {} surf;
        int cycles;
        if (p.z > sectors[s].z1) {
            sectors[s].surf = 1;
            cycles = 2;
            memset(sectors[s].slut, scrH, sizeof(int) * scrW);
        }
        else if (p.z < sectors[s].z2) {
            sectors[s].surf = 2;
            cycles = 2;
            memset(sectors[s].slut, 0, sizeof(int) * scrW);
        }
        else {
            sectors[s].surf = 0;
            cycles = 1;
        }


        for (int i = 0; i < cycles; i++) {
            for (w = sectors[s].idx; w < sectors[s].end; w++) {

            //     pixel(walls[w].x1 + 100, walls[w].y1 + 100, walls[w].c);
            //     pixel(walls[w].x2 + 100, walls[w].y2 + 100, walls[w].c);
                int wx[2],
                    wy[2],
                    wb[2],
                    wt[2];

                //translation
                int x1 = walls[w].x1 - p.x,
                    y1 = walls[w].y1 - p.y;
                
                int x2 = walls[w].x2 - p.x,
                    y2 = walls[w].y2- p.y;
                
                if (i) {
                    int swap;
                    swap = x1;
                    x1 = x2;
                    x2 = swap;

                    swap = y1;
                    y1 = y2;
                    y2 = swap;
                }

                //rotate
                wx[0] = x1 * p.ca - y1 * p.sa;
                wx[1] = x2 * p.ca - y2 * p.sa;

                wy[0] = y1 * p.ca + x1 * p.sa;
                wy[1] = y2 * p.ca + x2 * p.sa;
                

                wb[0] = sectors[s].z1 - p.z + p.l * wy[0];// / 32.0;
                wb[1] = sectors[s].z1 - p.z + p.l * wy[1];// / 32.0;


                // clip
                if (wy[0] < 1 && wy[1] < 1) {
                    continue;
                }
                if (wy[0] < 1) {
                    clipBehind(&wx[0], &wy[0], &wb[0], wx[1], wy[1], wb[1]);
                }
                if (wy[1] < 1) {
                    clipBehind(&wx[1], &wy[1], &wb[1], wx[0], wy[0], wb[0]);
                }

                // top z
                wt[0] = wb[0] + sectors[s].z2;
                wt[1] = wb[1] + sectors[s].z2;


                //screen space

                wx[0] = wx[0] * FOCAL / wy[0] + scrW/2;
                wx[1] = wx[1] * FOCAL / wy[1] + scrW/2;

                wt[0] = wt[0] * FOCAL / wy[0] + scrH/2;
                wt[1] = wt[1] * FOCAL / wy[1] + scrH/2;

                wb[0] = wb[0] * FOCAL / wy[0] + scrH/2;
                wb[1] = wb[1] * FOCAL / wy[1] + scrH/2;

                // if (wx[0] < 0 && wx[0] >= scrW) return;
                // if (wx[1] < 0 && wx[1] >= scrW) return;
                // if (wb[0] < 0 && wb[0] >= scrH) return;
                // if (wb[1] < 0 && wb[1] >= scrH) return;
                // if (wt[0] < 0 && wt[0] >= scrH) return;
                // if (wt[1] < 0 && wt[1] >= scrH) return;
                
                drawWall(wx[0], wx[1], wb[0], wb[1], wt[0], wt[1], s, w, i);
            }
            // sectors[s].surf *= -1;
        }
    }
}

int main() {
    loadMap("C:\\Users\\User\\Pictures\\doom-style-game\\map.txt");
    
	SDL_Init(SDL_INIT_VIDEO);

    SDL_Rect Bounds;
    SDL_GetDisplayBounds(0, &Bounds);
    scrW = 640;//Bounds.w;
    scrH = 480;//Bounds.h;

    window = SDL_CreateWindow("doom-style game!", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, scrW, scrH, SDL_WINDOW_SHOWN);
	surface = SDL_GetWindowSurface(window);
	pixels = (uint32_t*)surface->pixels;



    for (int s = 0; s < sectCount; s++) {
        //("%d %d %d %d", sectors[s].idx, sectors[s].end, sectors[s].z1, sectors[s].z2);
        //calculates the sector position as a vector
        double vx;
        double vy;
        for (int w = sectors[s].idx; w < sectors[s].end; w++) {
            vx += (walls[w].x2 + walls[w].x1)/2;
            vy += (walls[w].y2 + walls[w].y1)/2;
        }
        sectors[s].x = vx/(sectors[s].end - sectors[s].idx);
        sectors[s].y = vy/(sectors[s].end - sectors[s].idx);

        sectors[s].slut = malloc(sizeof(int) * scrW);
    }

    double deltaTime, lastFrameTime = SDL_GetTicks(), currentTime;
    #define P_SPD 0.2
    #define P_RSPD 0.002
 	bool quit = false;
	SDL_Event e;
	while (!quit) {
        currentTime = SDL_GetTicks();
        deltaTime = currentTime - lastFrameTime;
        lastFrameTime = currentTime;
        // printf("%f\n", deltaTime);
		while (SDL_PollEvent(&e)) {
			if (e.type == SDL_QUIT) {
				quit = true;
			}
			else if (e.type == SDL_KEYDOWN) {
				keys[e.key.keysym.scancode] = true;
			}
			else if (e.type == SDL_KEYUP) {
				keys[e.key.keysym.scancode] = false;
			}
		}

            p.sa = sin(p.a);
            p.ca = cos(p.a);
            //
            double movex = 0, movey = 0, movez = 0;
            if (keys[SDL_SCANCODE_W])
                movey++;
            if (keys[SDL_SCANCODE_S])
                movey--;
            if (keys[SDL_SCANCODE_A])
                movex--;
            if (keys[SDL_SCANCODE_D])
                movex++;
            if (keys[SDL_SCANCODE_SPACE])
                movez++;
            if (keys[SDL_SCANCODE_LSHIFT])
                movez--;
            if (keys[SDL_SCANCODE_LEFT])
                p.a -=  P_RSPD * deltaTime;
            if (keys[SDL_SCANCODE_RIGHT])
                p.a +=  P_RSPD * deltaTime;
            if (keys[SDL_SCANCODE_UP])
                p.l -=  P_RSPD * deltaTime;
            if (keys[SDL_SCANCODE_DOWN])
                p.l +=  P_RSPD * deltaTime;

                
            int oldx = movex;
            movex = movex * p.ca + movey * p.sa;

            movey = movey * p.ca - oldx * p.sa;

            p.x += movex * P_SPD * deltaTime;
            p.y += movey * P_SPD * deltaTime;
            p.z += movez * P_SPD * deltaTime;
        
		SDL_LockSurface(surface);
		memset(pixels, 0, sizeof(uint32_t) * scrW * scrH);
            render();
            pixel(p.x, p.y, 0xFFFFFFFF);
		SDL_UnlockSurface(surface);
		SDL_UpdateWindowSurface(window);
    }
    for (int s = 0; s < sectCount; s++) {
        free(sectors[s].slut);
    }
	SDL_DestroyWindow(window);
	SDL_Quit();
    return 0;
}