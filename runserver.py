import os
from biotin_flask import app

if __name__ == '__main__':
    # General Testing
    #app.run(debug=True)

    # Heroku Deployment
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)

    # Misc
    #app.run(host=os.getenv('IP', '0.0.0.0'), port=int(os.getenv('PORT', 8080)), debug=True)
    #app.run(host=os.getenv('IP', '0.0.0.0'), port=int(os.getenv('PORT', 8080)), debug=False)
    #app.run()
