#include <QtCore/QCoreApplication>
#include <Qt/qimage.h>
#include <Qt/qpainter.h>

#include <math.h>
#include <stdio.h>

#define MAX_FACTOR      0x40000000

typedef unsigned long  uint32;

uint32 prime_table[MAX_FACTOR / (sizeof(uint32) * 8) + 2];

int is_prime(uint32 num)
{
  if(num > MAX_FACTOR)
      return 1;
  return (prime_table[num / (sizeof(uint32) * 8)] &
                  (1L << (num % (sizeof(uint32) * 8))));
}

void unset_prime(uint32 num)
{
  if(num > MAX_FACTOR)
      return;
  prime_table[num / (sizeof(uint32) * 8)] &=
      ~(1L << (num % (sizeof(uint32) * 8)));
}

void set_prime(uint32 num)
{
  if(num > MAX_FACTOR)
      return ;
  prime_table[num / (sizeof(uint32) * 8)] |=
      (1L << (num % (sizeof(uint32) * 8)));
}

QPoint point(double theta) {
    theta += 2 * 4 * atan(1);

    int    block  = 0;
    double _local = theta;
    while(_local > 2 * 4 * atan(1)) {
        int step = log(_local) / log(2 * 4 * atan(1));
        block  += step;
        _local -= step * 2 * 4 * atan(1);
    }
    int    local  = (int)(_local / (2 * 4 * atan(1)) * 6 * block);
    local = local <= 3 * block ? local : local - 6 * block;

    return QPoint(block - abs(local),
                  abs(local) > 3 * block / 2 ?
                      local :
                      (local < 0 ? -1 : 1) * (3 * block / 2 - abs(local))
                      );
}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    // make prime.
    uint32 n_primes;
    unset_prime(0);
    unset_prime(1);
    for(uint32 i = 2; i < MAX_FACTOR; i ++)
        set_prime(i);
    n_primes = 0;
    for(uint32 i = 2; i < MAX_FACTOR; i ++)
        if(is_prime(i)) {
            for(uint32 j = i * 2; j < MAX_FACTOR; j += i)
                  unset_prime(j);
            n_primes ++;
        }

    // 2 min, 24fps, 10 sec / increment.
    const int frames_per_increment = 240000;
    const int frame_unit_start = 0;
    const int frame_unit_end   = 1;
    const int width  = 779;
    const int height = 779;

    if(a.arguments().count() < 2) {
        puts("ulam <directory> <start frame> <end frame>\n");
        puts("\t#frame number#.bmp will be generated.\n");
        return a.exec();
    }

    int i_start = frame_unit_start * frames_per_increment + 1;
    if(a.arguments().count() > 2)
        i_start = a.arguments().at(2).toInt() + frame_unit_start * frames_per_increment;

    int i_end = frame_unit_end * frames_per_increment;
    if(a.arguments().count() > 3)
        i_end = a.arguments().at(3).toInt() + frame_unit_start * frames_per_increment + 1;

    int frame_offset = 0;
    if(a.arguments().count() > 4)
        frame_offset = a.arguments().at(4).toInt();
  
    for(int i = i_start; i < i_end; i ++) {
        QImage img(QSize(width, height), QImage::Format_RGB32);
        QPainter p;

        int starts[width][height];
        int ends[width][height];

        for(int x = 0; x < width; x ++)
            for(int y = 0; y < height; y ++) {
                starts[x][y] = 0;
                ends[x][y] = 0;
            }

        uint32 max = 0;
        for(uint32 j = 0; j < MAX_FACTOR; j ++) {
            for(; j < MAX_FACTOR; j ++)
                if(is_prime(j)) break;

            double r = j * (double)i / (double)frames_per_increment;
            QPoint pt(point(r));

            if(pt.x() * pt.x() + pt.y() * pt.y() >
                    (width / 2) * (width / 2) + (height / 2) * (height / 2))
                break;

            int x = pt.x() + width / 2;
            int y = pt.y() + height / 2;

            if(x < 0 || width <= x || y < 0 || height <= y)
                continue;

            if(starts[x][y] == 0)
                starts[x][y] = j;
            ends[x][y] = j;

            max = j;
        }

        p.begin(&img);

        for(int x = 0; x < width; x ++)
            for(int y = 0; y < height; y ++) {
                QColor color;
                color.setHsl((starts[x][y] < 3 ? 0 : 255. * log(starts[x][y]) / log(max)),
                             (ends[x][y] < 3 ? 0 : 255. * log(ends[x][y]) / log(max)),
                             128);
                p.setPen(color);
                p.drawPoint(x, y);
            }
        p.end();

        img.save(a.arguments().at(1) + QString("/") +
                 QString("%1").arg((frame_unit_end - frame_unit_start) * frames_per_increment -
                                   (i - frame_unit_start * frames_per_increment) - frame_offset)
                 + QString(".bmp"), "bmp");
    }

    // something bugly here.
    //    return a.exec();
    return 0;
}
