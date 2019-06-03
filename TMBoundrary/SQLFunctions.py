import psycopg2
import psycopg2.extras


class RowSQL:
    _conn_param_dict_def = {
        "dbname": "ecod",
        "user": "ecodweb",
        "host": "129.112.32.63",
        "port": "45000",
        "password": "serveruse421",
    }
    def __init__(self, dbname: str = None, user: str = None,
                 host: str = None, port: str = None, password: str = None):
        self._conn_param_dict = dict()
        print(self._conn_param_dict_def.keys())
        if not dbname is  None:
            self._conn_param_dict["dbname"] = dbname
        else:
            self._conn_param_dict["dbname"] = self._conn_param_dict_def['dbname']
        if not user is  None:
            self._conn_param_dict["user"] = user
        else:
            self._conn_param_dict["user"] = self._conn_param_dict_def['user']
        if not host is  None:
            self._conn_param_dict["host"] = host
        else:
            self._conn_param_dict["host"] = self._conn_param_dict['host']
        if   not port is  None:
            self._conn_param_dict["port"] = port
        else:
            self._conn_param_dict["port"] = self._conn_param_dict['port']
        if not password is  None:
            self._conn_param_dict["password"] = password
        else:
            self._conn_param_dict["password"] = self._conn_param_dict['password']

    def _connect(self):
        self._conn = psycopg2.connect(**self._conn_param_dict)

    def get_domain_row(self,domain_id: str):
        cur = self._conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
        cur.execute('SELECT * FROM domain where id=%s', (domain_id,))
        results = cur.fetchone()
        if results is None:
            print(f"cannot find {domain_id}")
        res_dict = {k: v for k, v in results.items()}
        cur.close()
        return res_dict
