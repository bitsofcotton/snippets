#include <QtCore/QCoreApplication>
#include <Qt/qimage.h>
#include <Qt/qpainter.h>

#include <math.h>
#include <stdio.h>

const double max2 = 16.;
const int max_count = 64;

int seeConverge(double x0, double y0, double k) {
    int cnt;
    double x = x0, y = y0;
    for(cnt = 0; cnt < max_count; cnt ++) {
        double r2 = x * x + y * y;
        if(r2 > max2 || !isfinite(r2)) return cnt;

        double r = exp(k * log(sqrt(r2)));

        // also have interests in atan2(x, y);
        double theta = k * atan2(y, x);

        // this is z = z^k + z_0, but also have interests in z = z^k + z loop.
        x = r * cos(theta) + x0;
        y = r * sin(theta) + y0;
    }
    return cnt;
}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    // 2 min, 24fps, 10 sec / increment.
    const int frames_per_increment = 240;
    const int frame_unit_start = -6;
    const int frame_unit_end   = 6;
    const int width  = 779;
    const int height = 779;

    const double view_width  = 6.;
    const double view_height = 6.;

    const double inc_x_unit = view_width / width;
    const double inc_y_unit = view_height / height;

    if(a.arguments().count() < 2) {
        puts("mandel <directory> <start frame> <end frame>\n");
        puts("\t#frame number#.bmp will be generated.\n");
        return a.exec();
    }

    int i_start = frame_unit_start * frames_per_increment;
    if(a.arguments().count() > 2)
        i_start = a.arguments().at(2).toInt() + frame_unit_start * frames_per_increment;

    int i_end = frame_unit_end * frames_per_increment;
    if(a.arguments().count() > 3)
        i_end = a.arguments().at(3).toInt() + frame_unit_start * frames_per_increment + 1;

    for(int i = i_start; i < i_end; i ++) {
        QImage img(QSize(width, height), QImage::Format_RGB32);
        QPainter p;
        p.begin(&img);

        int    xidx, yidx;
        double x, y;
        for(xidx = 0, x = - view_width / 2;
            xidx < width;
            xidx ++, x += inc_x_unit) {
            for(yidx = 0, y = - view_height / 2;
                yidx < height;
                yidx ++, y += inc_y_unit) {
                // also have interest in z0 = 1/z, (x0, y0) = (x / (x^2 + y^2) - y / (x^2 + y^2)).
                int cnt = seeConverge(x, y, i / (double)frames_per_increment);

                QColor color;
                double r = cnt / (double)max_count;
                color.setHsl((r < 1. / 2. ? 2. * r : 1.) * 255,
                             (r >= 1. / 2. ? 2. * (1 - r) : 1.) * 255,
                              128);
                p.setPen(color);
                p.drawPoint(QPoint(xidx, yidx));
            }
        }
        p.end();

        img.save(a.arguments().at(1) + QString("/") +
                 QString("%1").arg(i - frame_unit_start * frames_per_increment)
                 + QString(".bmp"), "bmp");
    }

    // something bugly here.
    //    return a.exec();
    return 0;
}
