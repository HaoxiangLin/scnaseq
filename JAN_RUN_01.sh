WORKDIR=./example/work_dir
SCNASEQ_DIR=$(pwd)
JAR_FILE=$SCNASEQ_DIR/target/scnaseq0.5.jar
DATA_DIR=/Users/ALAN/Documents/git/scnaseq/data
INPUT=HCC1143_M101vN.counts
EXPERIMENT_ID=TEST_JAN

echo $WORKDIR
echo $JAR_FILE
cd $WORKDIR

#--------------------------------------
# MODES
#--------------------------------------
#sCNA-Seq is run with 3 different modes.
#MODE 1 = Whole genome training - Defines cellularity and ratio - uses large windows, low copy number states, high number of interations
#MODE 2 = Whole genome characterisation - uses cellularity and ratio from MODE 1 to characterise genome at higher resolution
#MODE 3 = Regional characterisation - able to provide ultra high resolution analysis of single region

#---------------------------------------
MODE=3
#--------------------------------------

#--------------------------------------
#Additional Mode Data
#--------------------------------------
#For Modes 2 and 3, you will also need to define the cellularity and the ratio
CELL="0.936"
RATIO="0.992"
COLON=":"
CellRatio=${CELL}${COLON}${RATIO}
echo "$CellRatio"

#-------------------------------------
# Regional Analysis
#-----------------------------------
#To run a regional analysis you will need to define the region you want to analyse
CHR_FROM="12"
START="0"
CHR_TO="12"
END="250"

RANGE="all:${CHR_FROM},${START}mb:${CHR_TO},${END}mb"
RANGE1="${CHR_FROM}""_${START}mb""_${CHR_TO}""_${END}mb"

#--------------------------------------
# MODE INTERALS
#--------------------------------------

if [ "$MODE" = "1" ]; then
        echo "WHOLE GENOME TRAINING MODE - qCNV window = $WINDOW - sCNAHitSeq Window = $cumR"
        cumR=1000
        MAXCN=2
        AVERAGE_DEPTH="false"
        MAKE_NEW_REF="false"
        train=2
        plot=1
        CellRatio="0.9999999:1.0"
        RANGE="all:1,0mb:X,300mb"
        RANGE1="1_0mb_X_300mb"
        BDW="10" #should be 10
        NumberIT=25
        echo "Determining Cellularity and CNA Ratio from tumour / normal genomes"
fi

if [ "$MODE" = "2" ]; then
        echo "WHOLE GENOME ANALYSIS MODE - sCNAHitSeq Window = $cumR - Cell/Ratio ${CellRatio}"
        cumR=100
        AVERAGE_DEPTH="true"
        MAKE_NEW_REF="false"
        train=0
        MAXCN=5
        plot=1
        RANGE="all:1,0mb:X,300mb"
        RANGE1="1_0mb_X_300mb"
        BDW="5" #should be 10
        NumberIT=1
        echo "Determining Cellularity and CNA Ratio from tumour / normal genomes"
fi

if [ "$MODE" = "3" ]; then
        cumR=2
        echo "REGIONAL sCNA CHARACTERISATION MODE - qCNV window = $WINDOW - sCNAHitSeq Window = $cumR"
        echo "User defined cellularity = $CELL, User defined sCNA ratio = $Ratio"
        AVERAGE_DEPTH="true"
        MAKE_NEW_REF="false"
        MAXCN=5
        train=0
        plot=1
        BDW="1"
        NumberIT=1
        echo "Determining sCNAs with user defined cellularity and ratio information, determined using mode 1 sCNAHit-Seq "

fi

#---------------------------------------------------
#Additional Variables
#--------------------------------------------------
ploidy=1
BC=2
#MAXCN=5
THREADS=1
CHROMS="1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22:X"
WINDOW=1000
BETA=2

java -cp  $JAR_FILE lc1.dp.appl.CNVHap --paramFile  param.txt --index 1 --column 1 --baseDir  $DATA_DIR   --numIt ${NumberIT}  --plot $plot  --showScatter false  --r_panel true  --b_panel true  --showHMM false    --showHaps false  --restrictKb 0kb:0kb  --include $INPUT  --experiment $EXPERIMENT_ID --initialCellularity ${CellRatio}  --downsample 1:1  --dilute 1:0   --cumulativeR $cumR   --useAvgDepth ${AVERAGE_DEPTH} --makeNewRef ${MAKE_NEW_REF}  --maxPloidy1 1  --noCopies ${ploidy}  --backgroundCount ${ploidy}   --ratioAsLevels null  --trainCellularity ${train}    --betaDownWeight ${BDW}   --noSampleFromBeta ${BETA}   --minNormalDepth 0 --mid ${RANGE}   --backgroundCount1 ${BC}  --numThreads 1 --alleles 1   --globalTrans false   --trainTransitions false --plotCNAsShape false --depthThresh 1  --onlyGlobalTrans  false --continueThresh 0 --reference NORMAL

