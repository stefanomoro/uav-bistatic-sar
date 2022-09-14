# DATASET

Files should be copied in folders with name composed like `YYYYMMDD`.
The complete folder structure should be like this

```
20220826
└───drone track
|   └───raw
│   │   flight1.mat
│   │   flight2.mat
└───radar
|   └───raw
|   |   |   test1.bb
|   |   |   test1_meta.mat
|   |   ...
|   └───RC
|   |   └───complete
|   |   |   |   test1.bb
|   |   |   |   test1_meta.mat
|   |   └───cut
|   |   |   |   test1_cut.bb
|   |   |   |   test1_cut_meta.mat
```

`RC` folder is the most important for the elaboration. In `complete` there are the file with all the 150km ranges, while in `cut` the file was cutted to only 1km range.
