#ls hg19*.conf|xargs -i echo "../src/redis-server {}"|bash
nohup ./src/redis-server aof_conf/hg19.0.conf & 
nohup ./src/redis-server aof_conf/hg19.10.conf &
nohup ./src/redis-server aof_conf/hg19.11.conf &
nohup ./src/redis-server aof_conf/hg19.1.conf & 
nohup ./src/redis-server aof_conf/hg19.2.conf & 
nohup ./src/redis-server aof_conf/hg19.3.conf &
nohup ./src/redis-server aof_conf/hg19.4.conf &
nohup ./src/redis-server aof_conf/hg19.5.conf &
nohup ./src/redis-server aof_conf/hg19.6.conf &
nohup ./src/redis-server aof_conf/hg19.7.conf &
nohup ./src/redis-server aof_conf/hg19.8.conf &
nohup ./src/redis-server aof_conf/hg19.9.conf &
