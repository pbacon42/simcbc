## Copyright (C) 2010 Eric Chassande-Mottin, CNRS (France)
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see
## <http://www.gnu.org/licenses/>.

## general routine to read ASCII data from a file and store into a cell struct array

function out=readcatalog(file,fields,delimiter,deblank)

  error(nargchk(2,4,nargin));
  
  if (nargin<3)
    delimiter=" ";
  endif

  if (nargin<4)
    deblank=false;
  endif

  if !exist(file,"file")
    error("readcatalog: file %s not found",file);
  endif

  descr=dir(file);
  if descr.bytes == 0
    warning("readcatalog: file %s is empty",file);
    out=cell;
    return
  endif
  
  ## remove comments
  if deblank
    deblank_cmd="| sed \"s/^ *//;s/ *$//;s/ \\{1,\\}/ /g\"";
  else
    deblank_cmd="";
  endif

  tmpfile=sprintf("/tmp/readcat%d.txt",unidrnd(1000));

  command=sprintf("grep -v \"^#\" %s %s > %s",file,deblank_cmd,tmpfile);
  system(command);

  array=csv2cell(tmpfile,delimiter);

  if size(array,2)!=length(fields)
    error("file %s has %d columns while %d fields are requested",file,size(array,2),length(fields));
  endif

  delete(tmpfile);

  out=cell2struct(array,fields,2);

endfunction
