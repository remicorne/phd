import os
import pandas as pd
import json
import sys
from io import StringIO

import plotly
import plotly.express as px

from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for
)
from werkzeug.exceptions import abort
from werkzeug.utils import secure_filename
from webApp.flaskr.auth import login_required

from module.quantitative import (
    justStats
)

from module.utils import subselectDf

from module.getters import (
    getQuantitativeStats,
)


ALLOWED_EXTENSIONS = {'csv', 'json'}
UPLOAD_FOLDER = os.path.join(os.getcwd(),'webApp','input')

bp = Blueprint('dashboard', __name__)


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@bp.route('/')
@login_required
def index():
   
    return render_template('blog/dashboard.html', data = {"rawfile":"","graphJSON":""})



@bp.route('/experimentalInfo', methods=['POST'])
@login_required
def upload_exp_info():

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


    
    return render_template('blog/dashboard.html', data = {"rawfile":"" ,"graphJSON":""})


@bp.route('/treatmentMapping', methods=['POST'])
@login_required
def upload_treatmentMapping():

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

    return render_template('blog/dashboard.html', data = {"rawfile":"" ,"graphJSON":""})



@bp.route('/data', methods=['POST'])
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


    data_canada = px.data.gapminder().query("country == 'Canada'")
    fig = px.bar(data_canada, x='year', y='pop')
    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

    justStats(filename, 
            experiments=['dose_response'], 
            compounds=['5HIAA/5HT'], 
            regions=["OF","PL","CC", "IC","M", "SJ","SL1", "SL6", "SR6", "SR1", "AC", "V",  
                    "Am", "dH", "vH", "NAc", "VM", "DM","VL", "DL", "MD",  "VPL",  "VPR", 
                    "DG", "Y",  "SC","SN", "VTA", "DR","MR", "CE"], 
            p_value_threshold=0.05)

    print("Calculation done!")

    # subselectDf(getQuantitativeStats(filename), {'experiment':'dose_response', 
    #                                             #  'is_significant':True, 
    #                                             'compound':'5HIAA/5HT', 
    #                                             'test':'one_way_anova',
    #                                             'region':["VPR"]}
    #         )

    print("Second step done !")

    return render_template('blog/dashboard.html', data = {"rawfile":"" ,"graphJSON":graphJSON})
    
   





    


    # with open(os.path.join(UPLOAD_FOLDER, filename),"r") as fid:
    #     data = fid.read()

    # df = pd.DataFrame({
    #     "Fruit": ["Apples", "Oranges", "Bananas", "Apples", "Oranges", "Bananas"],
    #     "Amount": [4, 1, 2, 2, 4, 5],
    #     "City": ["SF", "SF", "SF", "Montreal", "Montreal", "Montreal"]
    # })

    # fig = px.bar(df, x="Fruit", y="Amount", color="City", barmode="group")