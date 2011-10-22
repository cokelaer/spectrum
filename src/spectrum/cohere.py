
def cohere(kwargs):

    if len(kwargs<2):
        raise ValueError("cohere: Need at least 2 args. Use help cohere.")

    nvarargin = len(kwargs)
"""
  %% remove any pwelch RESULT args and add 'trans'
  for iarg=1:nvarargin
    arg = varargin{iarg};
    if ( ~isempty(arg) && ischar(arg) && ( strcmp(arg,'power') || ...
           strcmp(arg,'cross') || strcmp(arg,'trans') || ...
           strcmp(arg,'coher') || strcmp(arg,'ypower') ))
      varargin{iarg} = [];
    end
  end
  varargin{nvarargin+1} = 'coher';
  %%
  saved_compatib = pwelch('R11-');
  if ( nargout==0 )
    pwelch(varargin{:});
  elseif ( nargout==1 )
    Pxx = pwelch(varargin{:});
    varargout{1} = Pxx;
  elseif ( nargout>=2 )
    [Pxx,f] = pwelch(varargin{:});
    varargout{1} = Pxx;
    varargout{2} = f;
  end
  pwelch(saved_compatib);
  saved_compatib = 0;
end
"""
