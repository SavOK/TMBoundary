import psycopg2
import psycopg2.extras

conn_param_dict = {
        "dbname":"ecod", 
        "user":"ecodweb", 
        "host":"129.112.32.63", 
        "port":"45000", 
        "password":"serveruse421",
        }

class RowSQL:
    def __init__(self, dbname:str=None, user:str=None, host:str=None, 
        port:str=None,password:str=None ):
        self._conn_param_dict = dict()
        if dbname is not None:
            self._conn_param_dict["dbname"]=dbname
        if user is not None:
            self._conn_param_dict["user"]=user
        if host is not None:
            self._conn_param_dict["host"]=host
        if port is not None:
            self._conn_param_dict["port"]=port
        if password is not None:
            self._conn_param_dict["password"]=password
        }
    
    def _connect(self):
        self._conn =  psycopg2.connect(**self._conn_param_dict)

def get_domain_row(domain_id:str, c):
    c.execute('SELECT * FROM domain where id=%s', (domain_id,))
    results = c.fetchone()
    if results is None:
        print(f"cannot find {domain_id}")
    res_dict = {k:v for k,v in results.items()}
    return res_dict