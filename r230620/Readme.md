- versions:
    - bmibdcaller_simulations: 1922a9daf7e6041a5e32303524e708337f109ef3
    - hmmibd-rs: c47b07005c67b7cf0cc3f32ea6c0d4d1d86e3a3a
    - ibdutils: 26ef401388df7d299d7a2df290115be0b415ca0c
    - tskibd: c2f703b561c35233ded108af00df3b3588046704

- environment:
```
chib_conda
conda activate bmibdcaller_simulations
```

- command:
`nextflow ../main.nf -profile sge -resume`

