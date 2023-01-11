function [D] =fast_marching(mask)
% Compute an approximation of the distance D of pixels
% such that mask(x)=0 to the region mask=1


  [nbc nbl]=size(mask(:,:,1));
  D=ones(nbc,nbl)*1000000.0;

  for i=1:nbc,
  for j=1:nbl,
      if(mask(i,j)==1)
          D(i,j)=0;
      end
  end
  end


  for j=1:nbl,
     for i=1:nbc,



	  if(i>1)
	    if(D(i-1,j)+1.0<D(i,j))
	      D(i,j)=D(i-1,j)+1.0;
        end
      end


	  if(j>1)
	    if(D(i,j-1)+1.0<D(i,j))
	      D(i,j)=D(i,j-1)+1.0;
        end
      end


	  if(i>1 && j>1)
	    if(D(i-1,j-1)+sqrt(2.0)<D(i,j))
	      D(i,j)=D(i-1,j-1)+sqrt(2.0);
        end
      end


	  if(i<nbc && j>1)
	    if(D(i+1,j-1)+sqrt(2.0)<D(i,j))
	      D(i,j)=D(i+1,j-1)+sqrt(2.0);
        end
      end

     end
  end



  for i=nbc-1:-1:1 
	if(D(i+1,j)+1.0<D(i,j))
	  D(i,j)=D(i+1,j)+1.0;
    end
  end



  for j=nbl:-1:1,
    for i=nbc:-1:1, 

	  if(i<nbc)
	    if(D(i+1,j)+1.0<D(i,j))
	      D(i,j)=D(i+1,j)+1.0;
        end
      end

	  if(j<nbl)
	    if(D(i,j+1)+1.0<D(i,j))
	      D(i,j)=D(i,j+1)+1.0;
        end
      end


	  if(i<nbc && j<nbl)
	    if(D(i+1,j+1)+sqrt(2.0)<D(i,j))
	      D(i,j)=D(i+1,j+1)+sqrt(2.0);
        end
      end

	  if(i>1 && j<nbl)
	    if(D(i-1,j+1)+sqrt(2.0)<D(i,j))
	      D(i,j)=D(i-1,j+1)+sqrt(2.0);
        end
      end

	

    end
  end

      for i=2:nbc, 
        if(D(i-1,j)+1.0<D(i,j))
            D(i,j)=D(i-1,j)+1.0;
        end
      end







end
