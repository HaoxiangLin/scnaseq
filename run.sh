WORKDIR=./example/work_dir
SCNASEQ_DIR=$(pwd)
JAR_FILE=$SCNASEQ_DIR/target/scnaseq0.5.jar
DATA_DIR=/home/lachlan/Data
INPUT=HCC1143_M101vN.counts
EXPERIMENT_ID=TEST

echo $WORKDIR
echo $JAR_FILE
cd $WORKDIR

java -cp  $JAR_FILE lc1.dp.appl.CNVHap --paramFile  param.txt --index 1 --column 1 --baseDir  $DATA_DIR   --numIt  20  --plot 0  --showScatter false  --r_panel true  --b_panel true  --showHMM false    --showHaps false  --restrictKb 0kb:0kb  --include $INPUT  --experiment $EXPERIMENT_ID --initialCellularity 0.999:1.0  --downsample 1:1  --dilute 1:0   --cumulativeR 100   --useAvgDepth true     --makeNewRef false  --maxPloidy1 1  --noCopies 1  --backgroundCount 1   --ratioAsLevels null  --trainCellularity 2    --betaDownWeight 1   --noSampleFromBeta 2   --minNormalDepth 0 --mid all:1,0mb:21,300mb  --backgroundCount1  2  --numThreads 1 --alleles 1   --globalTrans false   --trainTransitions false --plotCNAsShape true --depthThresh 1  --onlyGlobalTrans  false --continueThresh 0 --reference NORMAL


