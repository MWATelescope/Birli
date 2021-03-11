# aoflagger::Strategy.Run Binding issue

When writing the Strategy.Run binding, I found that aoflagger would ignore existing flags passed in to it.

Here is some C++ code to demonstrate this

```c++
#include <aoflagger.h>
#include <iostream>

using namespace aoflagger;
using namespace std;

void printFlagMask(FlagMask mask)
{
        bool* buffer = mask.Buffer();
        int height = mask.Height();
        int width = mask.Width();
        int stride = mask.HorizontalStride();
        printf("height(chans)=%2d, width(steps)=%2d, stride=%2d @%08p \n", height, width, stride, &mask);
        for( int y=0; y<height; y++ ) {
                for( int x=0; x<stride; x++ ) {
                        printf("%2x", buffer[y * stride + x]);
                }
                cout << endl;
        }
}

size_t WIDTH = 5;
size_t HEIGHT = 6;

/**
 * using the default MWA_TELESCOPE flagging strategy, and an ImageSet with arbitrary width, height,
 * count and initialValue, perform a flagging run and print the resulting FlagMask
 *
 * if existingFlagMask is provided, then provide it to the flagging strategy.
 */
void flaggingDemo(AOFlagger* aoflagger, FlagMask* existingFlagMask = nullptr)
{
        size_t count = 4;
        float initialValue = 5;

        ImageSet imageSet = aoflagger->MakeImageSet(WIDTH, HEIGHT, count, initialValue);
        /* pollute one value in the imageSet */
        size_t imageStride = imageSet.HorizontalStride();
        float* imageBuffer = imageSet.ImageBuffer(2);
        imageBuffer[4 * imageStride + 3] = 999.0;

        string strategyFilePath = aoflagger->FindStrategyFile(TelescopeId::GENERIC_TELESCOPE, "minimal");

        Strategy strategy = aoflagger->LoadStrategyFile(strategyFilePath);
        FlagMask freshFlagMask = strategy.Run(imageSet, *existingFlagMask);
        printFlagMask(freshFlagMask);
}


int main(int argc, char* argv[])
{
        cout << "aoflagger library version: " << AOFlagger::GetVersionString() << endl;

        AOFlagger* aoflagger = new AOFlagger();

        cout << "run with no existing flag mask:" << endl;
        flaggingDemo(aoflagger);

        FlagMask existingFlagMask = aoflagger->MakeFlagMask(WIDTH, HEIGHT, false);
        /* uint8_t* existingMaskBuffer = reinterpret_cast<uint8_t *>(existingFlagMask.Buffer()); */
        bool* existingMaskBuffer = existingFlagMask.Buffer();
        size_t maskStride = existingFlagMask.HorizontalStride();

        /* populate the maskBuffer with arbitrary data */
        existingMaskBuffer[0] = true;
        existingMaskBuffer[2] = true;
        existingMaskBuffer[3] = true;

        cout << "existing flag mask:" << endl;
        printFlagMask(existingFlagMask);

        cout << "run with existing flag mask:" << endl;
        flaggingDemo(aoflagger, &existingFlagMask);

        return 0;
}
```

run with

```bash
g++ --std=c++11 -o aoflagger-test aoflagger-test.cc -laoflagger && chmod +x ./aoflagger-test && LD_LIBRARY_PATH=/usr/local/lib ./aoflagger-test
```

output:

```txt
aoflagger library version: 3.0
run with no existing flag mask:
height(chans)= 6, width(steps)= 5, stride= 8 @0x7ffc233199f0
 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0
 0 0 0 1 0 0 0 0
 0 0 0 0 0 0 0 0
existing flag mask:
height(chans)= 6, width(steps)= 5, stride= 8 @0x7ffc23319ab0
 1 0 1 1 0 0 0 0
 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0
run with existing flag mask:
height(chans)= 6, width(steps)= 5, stride= 8 @0x7ffc233199f0
 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0
 0 0 0 1 0 0 0 0
 0 0 0 0 0 0 0 0
 ```

I would expect the flag output out `run with existing flag mask` to be equal to `existing flag mask` `OR` `run with no existing flag mask`, so:

```txt
 1 0 1 1 0 0 0 0
 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0
 0 0 0 1 0 0 0 0
 0 0 0 0 0 0 0 0
```
