# 02/06/2022 Experimens

## Load tests

These notes concern the testing of various components to be used as load.

We carried out  three experiments:

1. Water, Vaseline, fish oil
2. Vaseline and fish oil
3. Fish oil only

In each case, the the 'adjustements' section of the scanner software was used to measure different frequency peaks. The results are below:

| Load                      | center frequency |
| ------------------------- | ---------------- |
| Water, Vaseline, fish oil | 41.980883        |
| vaseline, fish oil        | 41.980741        |
| fish oil only             | 41.980742        |

- the difference between the water peak and fat peak is 142 Hz. The [expected difference](https://mriquestions.com/f-w-chemical-shift.html) would be 220/1.5 = 146 Hz, so we are bang on.
- There is 1 Hz difference between fish oil and Vaseline; as such **Vaseline is an appropriate load for a phantom based on fish oil capsules**. 

## New slice tests

- The new slice worked well
- We tested GRE and TSE, both forward and reverse gradients. 
- From this, we can conclude that when the same bandwidth and FOV is used, distortion is the same for TSE and GRE images. 
- Our software failed to find some images in the GRE image; I have to investigate this.
