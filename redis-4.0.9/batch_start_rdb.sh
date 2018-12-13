#ls hg19*.conf|xargs -i echo "../src/redis-server {}"|bash
nohup taskset -c 0 ./src/redis-server rdb_conf/hg19.0.conf & 
nohup taskset -c 2 ./src/redis-server rdb_conf/hg19.10.conf &
nohup taskset -c 4 ./src/redis-server rdb_conf/hg19.11.conf &
nohup taskset -c 6 ./src/redis-server rdb_conf/hg19.1.conf & 
nohup taskset -c 8 ./src/redis-server rdb_conf/hg19.2.conf & 
nohup taskset -c 10 ./src/redis-server rdb_conf/hg19.3.conf &
nohup taskset -c 12 ./src/redis-server rdb_conf/hg19.4.conf &
nohup taskset -c 14 ./src/redis-server rdb_conf/hg19.5.conf &
nohup taskset -c 16 ./src/redis-server rdb_conf/hg19.6.conf &
nohup taskset -c 18 ./src/redis-server rdb_conf/hg19.7.conf &
nohup taskset -c 20 ./src/redis-server rdb_conf/hg19.8.conf &
nohup taskset -c 22 ./src/redis-server rdb_conf/hg19.9.conf &
