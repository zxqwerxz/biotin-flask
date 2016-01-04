import os
from biotin_flask import app

if __name__ == '__main__':
    #app.run(debug=True)
    #app.run(host=os.getenv('IP', '0.0.0.0'), port=int(os.getenv('PORT', 8080)), debug=True)
    app.run(host=os.getenv('IP', '0.0.0.0'), port=int(os.getenv('PORT', 8080)), debug=False)
    #app.run()