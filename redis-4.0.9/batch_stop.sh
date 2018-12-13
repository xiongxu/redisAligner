ps -u `whoami` -f |grep 'redis-server'|grep -v grep|awk '{print "kill -9 "$2}'|bash
