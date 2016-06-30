#mvn install:install-file \
#	-Dfile=./bzip2.jar \
#	-DgroupId=org.apache.tools \
#	-DartifactId=bzip2 \
#	-Dversion=1.0 \
#	-Dpackaging=jar \
#	-DgeneratePom=true

mvn install:install-file \
	-Dfile=./lib/hadoop-common-2.4.0.jar \
	-DgroupId=org.apache.hadoop \
	-DartifactId=hadoop-common \
	-Dversion=2.4.0-local \
	-Dpackaging=jar \
	-DgeneratePom=true

mvn install:install-file \
	-Dfile=./lib/hadoop-mapreduce-client-core-2.4.0.jar \
	-DgroupId=org.apache.hadoop \
	-DartifactId=hadoop-mapreduce-client-core \
	-Dversion=2.4.0-local \
	-Dpackaging=jar \
	-DgeneratePom=true
