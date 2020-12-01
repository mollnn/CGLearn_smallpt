#include <bits/stdc++.h>

#include "ColorSpace.hpp"

#include "Timer.hpp"
#include "Random.hpp"
#include "Vector.hpp"
#include "Ray.hpp"
#include "Material.hpp"
#include "Sphere.hpp"
#include "Scene.hpp"
#include "PathTracing.hpp"

double Tent(double x)
{
    return x < 1 ? sqrt(x) - 1 : 1 - sqrt(2 - x);
}

int main(int argc, char *argv[])
{
    srand(time(NULL));

    Timer timer_progress;
    Timer timer_rendertime;

    // 设定渲染参数
    int img_w = 512;                                                      // 图像宽度
    int img_h = 384;                                                      // 图像高度
    int subpixel_per_pixel_w = 2;                                         // 宽度方向每个像素拆分成的子像素个数
    int subpixel_per_pixel_h = 2;                                         // 高度方向每个像素拆分成的子像素个数
    int sample_per_subpixel = 4;                                          // 每个子像素的采样次数
    Ray camera(Vector(50, 52, 295.6), Vector(0, -0.042612, -1).Normal()); // 相机的位置与方向
    double film_size = 0.5;                                               // 投影平面线度（高度）
    Vector film_w = Vector(img_w * film_size / img_h);                    // 投影平面的宽度向量
    Vector film_h = (film_w % camera.direction).Normal() * film_size;     // 投影平面的高度向量
    Vector *img_buffer = new Vector[img_w * img_h];                       // 图像缓冲区
    double clip_near = 140;                                               // 最小渲染距离

    // 渲染，遍历每个像素，执行光线跟踪
    for (int img_y = 0; img_y < img_h; img_y++)
    {
        for (unsigned short img_x = 0; img_x < img_w; img_x++)
        {
            int pixel_id = (img_h - img_y - 1) * img_w + img_x; // 像素编号
            // 遍历每个子像素（每个像素被分成 2*2 个子像素）
            for (int subpixel_y = 0; subpixel_y < subpixel_per_pixel_h; subpixel_y++)
            {
                for (int subpixel_x = 0; subpixel_x < subpixel_per_pixel_w; subpixel_x++)
                {
                    Vector ans; // 该子像素的计算结果
                    // 对每个子像素进行若干次采样
                    for (int subsample_id = 0; subsample_id < sample_per_subpixel; subsample_id++)
                    {
                        double rand1 = 2 * RAND();
                        double rand2 = 2 * RAND();
                        double disp_x = Tent(rand1); // 采样相对子像素中心的 X 偏移
                        double disp_y = Tent(rand2); // 采样相对子像素中心的 Y 偏移
                        Vector ray_direction = film_w * (((subpixel_x + .5 + disp_x) / subpixel_per_pixel_w + img_x) / img_w - .5) +
                                               film_h * (((subpixel_y + .5 + disp_y) / subpixel_per_pixel_h + img_y) / img_h - .5) + camera.direction;
                        ans = ans + PathTracing(Ray(camera.origin + ray_direction * clip_near, ray_direction.Normal()), 0) * (1.0 / sample_per_subpixel);
                    }
                    // 将贡献加入图像缓冲区中
                    img_buffer[pixel_id] = img_buffer[pixel_id] + Vector(RadianceToBrightness(ans.x), RadianceToBrightness(ans.y), RadianceToBrightness(ans.z)) * (1.0 / subpixel_per_pixel_h / subpixel_per_pixel_w);
                }
            }

            // 渲染进度计算与报告
            timer_progress.End();
            if (timer_progress.GetTime() > 1.0)
            {
                std::cout << "Rendering Progress: " << std::setiosflags(std::ios::fixed)
                          << std::setprecision(2) << (double)((img_y)*100.0 / img_h + (img_x + 1) * 100.0 / img_w / img_h) << "%" << std::endl;
                timer_progress.Start();
            }
        }
    }

    timer_rendertime.End();
    std::cout << "============" << std::endl;
    std::cout << "Finish !" << std::endl;
    std::cout << "Rendering Time: " << timer_rendertime.GetTime() << " sec" << std::endl;

    FILE *f = fopen(argv[1], "w"); // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", img_w, img_h, 255);
    for (int i = 0; i < img_w * img_h; i++)
        fprintf(f, "%d\n%d\n%d\n", ColorFloat2Int(img_buffer[i].x), ColorFloat2Int(img_buffer[i].y), ColorFloat2Int(img_buffer[i].z));
}
