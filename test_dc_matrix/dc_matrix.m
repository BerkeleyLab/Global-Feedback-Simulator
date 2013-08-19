%%% Function to generate a generalized response matrix according to
%%% the BBFB settings defined in ../top_level/bbf_config.m.
%%% Each row of bbf.measurements represents a particular output of
%%% double compress (e.g. Energy deviation, bunch length deviation,
%%% etc.). First row of bbf.corrections represents amplitude,
%%% second represents phases. The value of each element represents
%%% the linac where the measurements/correction is taken
%%% from/applied to.
%%% IMPORTANT: if one wants to expand the measurements set by
%%% including new outputs of double compress, meaning expanding the
%%% rows of M, they should pay extreme attention to expanding both
%%% the settings in ../top_level/bbf_config.m and the vectors
%%% storing double compress outputs (do not forget to also modify
%%% dc_out!)
%%%
%%% Jack Olivieri - LBNL/Februrary 2013

function [M] = dc_matrix(paramsdc, Q_nom, bbf)
    
  row = size(bbf.measurements,2);
  column = size(bbf.corrections,2);
  M = zeros(row,column);
  ddE = zeros(5,2);
  sz = zeros(5,2);
  dt = zeros(5,2);
  sd = zeros(5,2);

  for m=1:row
    m
      linac = bbf.measurements(bbf.idx_meas(m),m);
      [indx1,indx2] = find((bbf.corrections<=linac)&(bbf.corrections>0));
      indx2
      for k=1:length(indx2)
	k
          excitation=zeros(2,5,2);
          col_exc = bbf.corrections(indx1(k),indx2(k));
          if (indx1(k)==2)
              row_exc = 1;
          else
              row_exc = 2;
          end             
          excitation(row_exc,col_exc,:)=(row_exc==2)*bbf.dv + (row_exc==1)*bbf.dphi
          for c=[1:2]
              [Ipk,sz(:,c),ddE(:,c),sd(:,c),dt(:,c),sdsgn,kk,Eloss] = ...
                  double_compressxv(paramsdc,0,0,0,0,0,0,excitation(1,:,c),excitation(2,:,c),Q_nom);
          end
          dc_out = [(ddE(linac,2)-ddE(linac,1)); (sz(linac,2)-sz(linac,1)); (dt(linac,2)-dt(linac,1));  (sd(linac,2)-sd(linac,1))];
          M(m,indx2(k)) = dc_out(bbf.idx_meas(m))/(excitation(row_exc,col_exc,2) - excitation(row_exc,col_exc,1));
      end
      clear indx1 indx2;
      
end
