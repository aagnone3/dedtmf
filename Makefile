DST=dedtmf
#TAR=${DST}.tgz
ZIP=${DST}.zip

#MATLAB=/usr/bin/Matlab 
MATLAB=/Applications/MATLAB_R2010b.app/bin/matlab
# MacOS 64 bit (dpwe-macbook)
DEPLOYTOOL=/Applications/MATLAB_R2010b.app/bin/deploytool
# Linux 64 bit (hog)
#DEPLOYTOOL=/usr/local/MATLAB/R2010b/bin/deploytool 
# Linux 32 bit (cherry)
#DEPLOYTOOL=${MATLAB} -r deploytool

DEMOFILE=demo_${DST}

THUMB=${DST}_thumb.png

SRCS=${DEMOFILE}.m ${DST}.m frame.m ola.m rmsenv.m dedtmf.py

DATA=${THUMB} phonexample.wav out-ref.wav

DEMOHTML=html/${DEMOFILE}.html
DEMOINDEX=html/index.html

DSTINDEX=${DST}/index.html

all: dist

# Test the python version against an earlier version
test:
	./dedtmf.py phonexample.wav out.wav
	sndcmp out-ref.wav out.wav

${DEMOHTML}: ${DEMOFILE}.m ${SRCS} ${DATA} 
	${MATLAB} -r "spublish ${DEMOFILE}; exit"

# spublish publishes with sound files

${DEMOINDEX}: ${DEMOHTML}
	sed -e 's@<div class="content">@<a href="http://www.ee.columbia.edu/~dpwe/">Dan Ellis</a> : <a href="http://www.ee.columbia.edu/~dpwe/resources/">Resources</a>: <a href="http://www.ee.columbia.edu/~dpwe/resources/matlab/">Matlab</a>: <div class="content"> <IMG SRC="'${THUMB}'" ALIGN="LEFT" HSPACE="10">@' -e 's/amp;auml;/auml;/g' < ${DEMOHTML} > ${DEMOINDEX}

MCR_NAME=${DST}_${ARCH}
PRJ_NAME=${DST}_prj

dist: ${SRCS} ${DATA} ${DEMOINDEX}
	rm -rf ${DST}
	mkdir ${DST}
	cp -p html/* ${DST}
	rm ${DST}/${DEMOFILE}.html
	cp -p ${SRCS} ${DATA} ${DST}
	rm -f ${DST}/*~
	-rm-extended-attribs.sh ${DST}
	zip -r ${ZIP} ${DST}
	cp -p ${ZIP} ${DST}
	rsync -avz ${DST} hog.ee.columbia.edu:public_html/LabROSA/projects/
	rsync -avz ${DST} labrosa.ee.columbia.edu:/var/www/LabROSA/projects/
	rm -rf ${DST}
