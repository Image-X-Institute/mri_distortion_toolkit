ok.

So so far we know:

- for some of the images we ahve to use `Gz=Gz*-1`
- We are not really getting consistent results with our ability to correct images.

### Does an 'end to end' script allow all images to correct 'themselves?'

- Transverse: yes, except with `Gz=Gz*-1`
- All others work

###  OK. Do the harmonics we extract and save also work to correct themselves?

- Yes they all work.

### OK. So what harmonics should be the same?

| Scan orientation | Slice direction | Phase direction | Frequency direction |
| ---------------- | --------------- | --------------- | ------------------- |
| Transverse       | z               | **y**           | **x**               |
| Coronal          | y               | **x**           | **z**               |
| Sagital          | x               | **z**           | **y**               |

## OK. start again. Comparing 'rot' harmonics versus sag harmonics.

**Rot harmonics**

```
x_mean:  0.4 mm +-  0.3 mm. max  1.4 mm
y_mean:  0.2 mm +-  0.2 mm. max  1.8 mm
z_mean:  0.4 mm +-  0.3 mm. max  2.3 mm
```

**Sag harmonics**

```
x_mean:  0.5 mm +-  0.4 mm. max  1.9 mm
y_mean:  1.0 mm +-  0.7 mm. max  3.6 mm
z_mean:  0.7 mm +-  0.6 mm. max  3.2 mm
```

Y in particular is a lot worse. Z also isn't great.

**Y Harmonics:**

![](/run/user/1000/doc/95456213/G_y_sag_vs_rot.png)

So, one of the main differences here is that series two (sagital) has a relatively big A_1_0 term which is lacking in the rotated harmonics. There's also a smaller difference in A_1_1. 

So what is this A_1_0 term... It could be some sort of B0 artefact, but we **should** be removing those...

Anway let's see what happens when we manually remove it...

```
x_mean:  0.5 mm +-  0.3 mm. max  1.7 mm
y_mean:  0.9 mm +-  0.6 mm. max  2.8 mm
z_mean:  0.5 mm +-  0.5 mm. max  3.2 mm
```

That's better, particularly regarding the maximum distortion. It's also improved the z distortion quite a lot. what happens if we "fix" the A_1_1 term...

```
x_mean:  0.5 mm +-  0.3 mm. max  1.7 mm
y_mean:  0.7 mm +-  0.5 mm. max  2.4 mm
z_mean:  0.5 mm +-  0.5 mm. max  3.3 mm
```

**Gz** looks very similar.

**Gx** looks like this...

![](/run/user/1000/doc/94b3af8d/G_x_sag_v_rot.png)

However - note that this is the slice select direction, and we probably have B0 effects in here which we are not correcting for. Particularly I suspect the A_1_0 and B_1_1 terms which are reversed are due to this...

Having said that, if we do "fix" these two terms, we get this result:

```
x_mean:  0.3 mm +-  0.3 mm. max  1.6 mm
y_mean:  0.3 mm +-  0.4 mm. max  2.5 mm
z_mean:  0.4 mm +-  0.4 mm. max  2.0 mm
```

So. 

- some of the differences appear to be B0 effects in the slice encode direction. These should be eliminated in the 'average' gradients anyway...
- The other major differences were in y. In principle, the B0 effects should be removed in this direction, but it doesn't mean it's working perfectly

## Flip it around.

Ok, how do the 'rotated' harmonics go at correcting the 'sagital' images?

```
x_mean:  1.1 mm +-  0.9 mm. max  5.8 mm
y_mean:  1.5 mm +-  1.0 mm. max  4.3 mm
z_mean:  0.4 mm +-  0.4 mm. max  3.4 mm
```

Dreadful!

So. We have two identical scans, and the harmonics derived from each scan do a great job correcting **itself** but not a very good job correcting other images...



## Correcting of the original 6 images:

- Sag: do pretty well, except in the X (slice direction)
- Tra: do well, except in Z (slice direction)
- Cor: do very well except for Y
- These results are what we expect. So why the hell are the rotated results so poor.
  - Could the reconstruction be backwards somehow?


| Slice orientation | x_mean (max)  | y_mean (max)  | z_mean (max)  |
| ----------------- | ------------- | ------------- | ------------- |
| Transverse        | 0.5 (4.9)     | 0.3 (2.3)     | **0.7 (4.0)** |
| Coronal           | 0.4 (1.6)     | **0.8 (6.7)** | 0.4 (3.5)     |
| Sagital           | **0.9 (5.7)** | 0.2 (1.2)     | 0.3 (3.7)     |

I have bolded the slice direction in each case. We don't correct very well in this direction, because we don't take into account the effects of B0 on the slice gradient.

then when we switch to the rotated:

```
x_mean:  0.9 mm +-  0.7 mm. max  4.0 mm
y_mean:  0.9 mm +-  0.7 mm. max  3.2 mm
z_mean:  0.8 mm +-  0.7 mm. max  3.6 mm
```

We would expect X to be a bit shit. But not the others...

Ok. So we are sort of back at square one... what is so different about these harmonics.

```
mean distortion:  1.8 mm, std:  0.9 mm, Max:  4.8 mm
x_mean:  1.0 mm +-  0.8 mm. max  4.1 mm
y_mean:  1.0 mm +-  0.7 mm. max  3.4 mm
z_mean:  0.7 mm +-  0.6 mm. max  3.3 mm
```

Here's the results of using just the sagital:

```
mean distortion:  1.6 mm, std:  0.7 mm, Max:  3.7 mm
x_mean:  0.6 mm +-  0.4 mm. max  1.9 mm
y_mean:  1.1 mm +-  0.8 mm. max  3.4 mm
z_mean:  0.7 mm +-  0.6 mm. max  3.2 mm
```

for the record here are the unocrrected results:

```
mean distortion:  5.5 mm, std:  3.1 mm, Max:  14.7 mm
x_mean:  2.5 mm +-  2.4 mm. max  11.7 mm
y_mean:  2.2 mm +-  1.8 mm. max  8.1 mm
z_mean:  3.3 mm +-  2.4 mm. max  11.0 mm
```





it just seems like this really is the right result... what else could effect it?

- Markers moved between scans? 
  - ... actually, that would do it. But it's hard to believe that that's the case...
- gradient performance different for different scans?
  - then the sagital scan should still work...
- 
