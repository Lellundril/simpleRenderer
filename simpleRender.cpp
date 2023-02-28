#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <cstring>
#include <sstream>
#include "tgaimage.h"

constexpr int WIDTH = 512;
constexpr int HEIGHT = 512;

std::vector<std::string> split (std::string s, std::string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}

// points est un vector de int : chaque triplet représente les coordonnées x y z du sommet lu
// face est un vector de int: par triplet, chaque nombre étant l'indice d'un des sommet du triangle courant. 
// donc indice du sommet (x1,y1,z1) puis (x2,y2,z2) puis (x3,y3,z3) et rebelote pour le triangle suivant
void readObj(std::string filename, std::vector<float> &vertex, std::vector<float> &vertexToScale, std::vector<int> &faces) {
    std::ifstream file(filename);
    if (file.is_open())
    {
        std::string line;
        std::string word;
        std::vector<std::string> vec;
        float px, py, pz; // point coord
        int f1, f2, f3; // face index
        while (getline(file, line))
        {            
            if(line.find("v ") != std::string::npos){ // Si la ligne contient "v " (= est un vertex)
                std::stringstream ss(line);
                getline(ss, word, ' '); // Ignore first word ("v")

                getline(ss, word, ' '); // x coordinate
                px = std::stof(word);

                getline(ss, word, ' '); // y coordinate
                py = std::stof(word);

                getline(ss, word, ' '); // z coordinate
                pz = std::stof(word);

                vertex.push_back(px); // Put in the coordinates of the corresponding vector
                vertex.push_back(py);
                vertex.push_back(pz);
                vertexToScale.push_back(px); // Put in the coordinates of the corresponding vector
                vertexToScale.push_back(py);
                vertexToScale.push_back(pz);
            }

            if(line.find("f ") != std::string::npos){ // Si la ligne contient "f " (= est un triangle/face)
                std::stringstream ss(line);

                getline(ss, word, ' '); // Ignore first word ("f")

                getline(ss, word, ' '); // First point
                vec = split(word, "/");
                f1 = std::stoi(vec.at(0));

                getline(ss, word, ' '); // Second point
                vec = split(word, "/");
                f2 = std::stoi(vec.at(0));
                
                getline(ss, word, ' '); // Third point
                vec = split(word, "/");
                f3 = std::stoi(vec.at(0));
                
                // std::cout << "faces = " << f1 << ", " << f2 << ", " << f3 << std::endl;
                faces.push_back(f1);
                faces.push_back(f2);
                faces.push_back(f3);
            }
        }
    }
}

// Bresenham
void drawLine(int x0, int x1, int y0, int y1, TGAImage &image, TGAColor color){
    bool steep = false;
    if (std::abs(x0-x1)<std::abs(y0-y1)) {
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    if (x0>x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    int dx = x1 - x0;
    int dy = y1 - y0;
    int y = y0;
    int derror3 = 2 * std::abs(dy); // Quantité d'erreur par rapport au segment à chaque "pas" sur l'axe x
    // L'erreur que j'introduis dans mon sys si je n'incrémente pas y
    int error3 = 0;
    for (int x = x0; x <= x1; x++) {
        if (steep)
            image.set(y, x, color);
        else
            image.set(x, y, color);
        error3 += derror3;
        if (error3 > .5 * dx) {
            y += dy>0 ? 1 : -1;
            error3 -= 2*dx;
        }
    }
}

// Draws a rectangle from two given opposite corners
void drawRectangle (int x1, int y1, int x2, int y2, TGAImage &image, TGAColor color){
    drawLine(x1, x1, y1, y2, image, color);
    drawLine(x1, x2, y1, y1, image, color);
    drawLine(x1, x2, y2, y2, image, color);
    drawLine(x2, x2, y1, y2, image, color);
}

// Draws every pixel in "vertex"
void displayPointsCloud(std::vector<float> &vertex, TGAImage &framebuffer, TGAColor color){
    for(int i = 0; i < vertex.size(); i = i+3){
        framebuffer.set(vertex[i], vertex[i+1], color);
    }
}

// Draws every triangle
void displayTriangles(std::vector<float> &vertex, std::vector<int> &faces, TGAImage &framebuffer, TGAColor color){
    int index, x1, y1, x2, y2, x3, y3;
    for(int i = 0; i < faces.size(); i = i+3){
        index = (faces.at(i)-1)*3;
        x1 = vertex.at(index);
        y1 = vertex.at(index + 1);

        index = (faces.at(i+1)-1)*3;
        x2 = vertex.at(index);
        y2 = vertex.at(index + 1);

        index = (faces.at(i+2)-1)*3;
        x3 = vertex.at(index);
        y3 = vertex.at(index + 1);

        drawLine(x1, x2, y1, y2, framebuffer, color);
        drawLine(x2, x3, y2, y3, framebuffer, color);
        drawLine(x3, x1, y3, y1, framebuffer, color);
    }
}

// Gets the coordinates for each triangle (line 'f x/y/z' in .obj) and puts them in "triangles", in the same order as the obj file
void getTriangles(std::vector<float> &vertex, std::vector<int> &faces, std::vector<float> &triangles){
    int index;
    float x1, y1, z1, x2, y2, z2, x3, y3, z3;
    for(int i = 0; i < faces.size(); i = i+3){
        index = (faces.at(i)-1)*3;
        x1 = vertex.at(index);
        y1 = vertex.at(index + 1);
        z1 = vertex.at(index + 2);

        index = (faces.at(i+1)-1)*3;
        x2 = vertex.at(index);
        y2 = vertex.at(index + 1);
        z2 = vertex.at(index + 2);

        index = (faces.at(i+2)-1)*3;
        x3 = vertex.at(index);
        y3 = vertex.at(index + 1);
        z3 = vertex.at(index + 2);

        triangles.push_back(x1);
        triangles.push_back(y1);
        triangles.push_back(z1);
        triangles.push_back(x2);
        triangles.push_back(y2);
        triangles.push_back(z2);
        triangles.push_back(x3);
        triangles.push_back(y3);
        triangles.push_back(z3);
    }
}

// Checks the sign of x1y1 relative to the segment (p2,p3). Tells which "side" of the line p1 is
float sign (int x1, int y1, int x2, int y2, int x3, int y3)
{
    return (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);
}

// Checks if point p is in triangle p1p2p3 or not
bool pointInTriangle (int px, int py, int x1, int y1, int x2, int y2, int x3, int y3)
{
    float d1, d2, d3;
    bool has_neg, has_pos;

    d1 = sign(px, py, x1, y1, x2, y2);
    d2 = sign(px, py, x2, y2, x3, y3);
    d3 = sign(px, py, x3, y3, x1, y1);

    has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(has_neg && has_pos);
}

// Function to display the contents of a given vector (!LOOPS)
void displayVector(auto vec){
    for (const auto& i : vec){
        std::cout << i << std::endl;
    }
}

int main(int argc, char const *argv[])
{
    // Constants needed throughout the program
    const TGAColor white = {255, 255, 255, 255};
    const TGAColor green = {0, 255, 0, 255};
    const TGAColor red = {255, 0, 0, 255};
    const TGAColor blue = {0, 0, 255, 255};
    TGAImage framebuffer(WIDTH, HEIGHT, TGAImage::RGB);

    // Data initialization
    std::vector<float> vertex;
    std::vector<float> vertexToScale;
    std::vector<int> faces;
    std::vector<float> triangles;
    std::vector<float> trianglesToScale;
    readObj("../african_head.obj", vertex, vertexToScale, faces);
    
    std::vector<float> light {0.,0.,1.}; // Vecteur lumière
    float scalar; // Le produit scalaire entre le vecteur lumière et la normale
    TGAColor greyLevel = {255, 255, 255, 255};

    // Put pixel coordinates to scale with window size
    for (int i = 0; i < vertexToScale.size(); i = i+3){
        vertexToScale[i] = (vertexToScale[i] * WIDTH) / 2 + WIDTH/2;
        vertexToScale[i+1] = (vertexToScale[i+1] * HEIGHT) / 2 + HEIGHT/2;
    }

    // Display points cloud :
    //displayPointsCloud(vertexToScale, framebuffer, white);

    // Display triangles
    //displayTriangles(vertexToScale, faces, framebuffer, white);

    // Get all triangles in a vector
    getTriangles(vertex, faces, triangles);
    getTriangles(vertexToScale, faces, trianglesToScale);

    for (int i = 0; i < trianglesToScale.size(); i = i+9){
        float x1, y1, z1, x2, y2, z2, x3, y3, z3;
        int ix1, iy1, iz1, ix2, iy2, iz2, ix3, iy3, iz3, cx1, cy1, cx2, cy2;

        // Prendre les coordonnées du triangle courant
        x1 = triangles[i];
        y1 = triangles[i+1];
        z1 = triangles[i+2];
        x2 = triangles[i+3];
        y2 = triangles[i+4];
        z2 = triangles[i+5];
        x3 = triangles[i+6];
        y3 = triangles[i+7];
        z3 = triangles[i+8];

        // Avec ces coordonnées, calculer la normale du triangle
        // Vecteur normal du triangle courant
        std::vector<float> norm = {(y2-y1)*(z3-z1) - (z2-z1)*(y3-y1), (z2-z1)*(x3-x1) - (x2-x1)*(z3-z1), (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)};
        // Longueur du vecteur normal
        float len = std::sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
        for (int i : {0,1,2}) norm[i] /= len; // Diviser la normale par sa longueur (normalise le vecteur)

        // Calculer le produit scalaire entre la normale et le vecteur "lumière"
        scalar = (norm[0] * light[0]) + (norm[1] * light[1]) + (norm[2] * light[2]);
        // Produit entre le scalaire obtenu et l'opacité du blanc pour déterminer le niveau de gris à appliquer
        if (scalar < 0) continue; // Si le scalaire est négatif (triangle tourné vers l'opposé de la caméra), ne rien faire
        greyLevel = {scalar*255, scalar*255, scalar*255, 255};

        int *zbuffer = new int[WIDTH*HEIGHT];

        // Prendre les coordonnées à l'échelle du triangle courant
        ix1 = trianglesToScale[i];
        iy1 = trianglesToScale[i+1];
        iz1 = trianglesToScale[i+2];
        ix2 = trianglesToScale[i+3];
        iy2 = trianglesToScale[i+4];
        iz2 = trianglesToScale[i+5];
        ix3 = trianglesToScale[i+6];
        iy3 = trianglesToScale[i+7];
        iz3 = trianglesToScale[i+8];

        // Get the coordinate of the square encapsulating the triangle
        cx1 = std::min({ix1, ix2, ix3});
        cx2 = std::max({ix1, ix2, ix3});
        cy1 = std::min({iy1, iy2, iy3});
        cy2 = std::max({iy1, iy2, iy3});

        // Draws every encapsulating rectangle in red
        //drawRectangle(cx1, cy1, cx2, cy2, framebuffer, red);

        // Goes through every point in the encapsulating square and colors it if it's inside the triangle
        int iz;
        for (int iy = cy1; iy <= cy2; iy++){
            for (int ix = cx1; ix <= cx2; ix++){
                if (pointInTriangle(ix, iy, ix1, iy1, ix2, iy2, ix3, iy3)){
                    iz = 0;
                    iz += iz1 /* * barycentric x*/;
                    iz += iz2 /* * barycentric y*/;
                    iz += iz3 /* * barycentric z*/;
                    if (zbuffer[ix + iy*WIDTH] < iz){
                        zbuffer[ix + iy*WIDTH] = iz;
                    }
                    framebuffer.set(ix, iy, greyLevel);
                }
            }
        }

    }

    framebuffer.write_tga_file("framebuffer.tga");
    return 0;
}
