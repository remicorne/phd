import os
import pandas as pd
import json


import plotly
import plotly.express as px

from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for
)
from werkzeug.exceptions import abort
from werkzeug.utils import secure_filename
from webApp.flaskr.auth import login_required


ALLOWED_EXTENSIONS = {'txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif'}
UPLOAD_FOLDER = os.path.join(os.getcwd(),'webApp','input')

bp = Blueprint('dashboard', __name__)


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@bp.route('/')
@login_required
def index():
   
    return render_template('blog/dashboard.html', data = {"rawfile":"","graphJSON":""})


@bp.route('/', methods=['POST'])
@login_required
def upload_file():

    # check if the post request has the file part
    if 'file' not in request.files:
        flash('No file part')
        return redirect(request.url)
    file = request.files['file']
    # If the user does not select a file, the browser submits an
    # empty file without a filename.
    if file.filename == '':
        flash('No selected file')
        return redirect(request.url)
    filename = secure_filename(file.filename)
    file.save(os.path.join(UPLOAD_FOLDER, filename))
    with open(os.path.join(UPLOAD_FOLDER, filename),"r") as fid:
        data = fid.read()

    # df = pd.DataFrame({
    #     "Fruit": ["Apples", "Oranges", "Bananas", "Apples", "Oranges", "Bananas"],
    #     "Amount": [4, 1, 2, 2, 4, 5],
    #     "City": ["SF", "SF", "SF", "Montreal", "Montreal", "Montreal"]
    # })

    # fig = px.bar(df, x="Fruit", y="Amount", color="City", barmode="group")
        
    data_canada = px.data.gapminder().query("country == 'Canada'")
    fig = px.bar(data_canada, x='year', y='pop')


    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)



    return render_template('blog/dashboard.html', data = {"rawfile":data,"graphJSON":graphJSON})
    
   


