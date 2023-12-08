from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for
)
from werkzeug.exceptions import abort

from webApp.flaskr.auth import login_required


bp = Blueprint('dashboard', __name__)


@bp.route('/')
@login_required
def index():
    return render_template('blog/dashboard.html', user=g.user["username"])

