import sys 
from datetime import datetime
timestamp = datetime.now().strftime("%Y/%m/%d_%H:%M:%S").replace("/", "_").replace(":", "_")
print(sys.argv[1], sys.argv[2], sys.argv[3],sys.argv[4])
print(timestamp)