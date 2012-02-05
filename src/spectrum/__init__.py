#: version
__version__ = '0.4.5'
#: default number os samples used to compute FFT
default_NPSD = 4096

import arma
import burg
import cholesky
import correlation
import correlog
import covar
import criteria
##try:
#    import dpss
#    from dpss import *
#except:
#    pass
import eigen
import eigenfre
#import fastrls
import levinson
#import linear_prediction
import linalg
#import lms
import lpc
import minvar
import modcovar
import periodogram
import psd
import datasets
import toeplitz
import tools
import window
import waveform
import yulewalker
import transfer


from transfer import *
from arma import *
from burg import *
from cholesky import *
from correlation import *
from correlog import *
from covar import *
from criteria import *
from eigen import *
from eigenfre import *
from levinson import *
#from linear_prediction import *
from linalg import *
#from lms import *
from lpc import *
from minvar import *
from modcovar import *
from periodogram import *
from psd import *
from datasets import *
from tools import *
from toeplitz import *
from window import *
from waveform import *
from yulewalker import *

