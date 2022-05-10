function globdat = postProcessStep(dir_output,props,globdat)

% Carries out all post-processing after a step has converged

if rem(globdat.ts.step,props.postProc.onScreen) == 0 || ...
       rem(globdat.ts.step,props.postProc.printVTK) == 0 || ...
       globdat.active == 0
       
    geometry       = props.geom;
    material       = props.mat;
    step           = globdat.ts.step;
    lodi           = globdat.lodi;
    tdisp          = globdat.state;
    PostProc       = props.postProc;
    PHTelem        = globdat.PHTelem;
    sizeBasis      = globdat.mesh.sizeBasis;
    numberElements = globdat.mesh.numElements;
    controlPts     = globdat.ctrlPts;
    errvec         = globdat.errVec{globdat.ts.step};

    dim            = geometry.dim;
    p              = geometry.p;
    q              = geometry.q;
    vIEN           = zeros(numberElements,4);
    physcoord      = zeros(4*numberElements,2);
    dispcoord      = zeros(4*numberElements,2);
    wcoord         = zeros(4*numberElements,1);
    sigmacoord     = zeros(4*numberElements,1);

    fudge = 0;
    nx = 2;
    ny = 2;
    px = linspace(-1+fudge,1-fudge, nx);
    py = linspace(-1+fudge,1-fudge, ny);

    %1D bernstein polynomials evaluated at the Gauss points on the master element
    [B_u,dB_u] = bernstein_basis(px,p);
    [B_v,dB_v] = bernstein_basis(py,q);

    dBdu = zeros(nx, ny, (p+1)*(q+1));
    dBdv = zeros(nx, ny, (p+1)*(q+1));
    R    = zeros(nx, ny, (p+1)*(q+1));

    % The derivatives of the 2D Bernstein polynomials at Gauss points on the master element
    basisCounter = 0;
    for j=1:q+1
        for i=1:p+1
            basisCounter = basisCounter + 1;
            dBdu(:,:,basisCounter) = dB_u(:,i)*B_v(:,j)';
            dBdv(:,:,basisCounter) = B_u(:,i)*dB_v(:,j)';
            R(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
        end
    end


    elementCounter = 0;
    for indexPatch = 1:length(PHTelem)
        for i=1:length(PHTelem{indexPatch})
            if isempty(PHTelem{indexPatch}(i).children)
                elementCounter =  elementCounter+1;

                vIEN(elementCounter,:) = [(elementCounter-1)*4+1:(elementCounter-1)*4+4];
                ximin = PHTelem{indexPatch}(i).vertex(1);
                ximax = PHTelem{indexPatch}(i).vertex(3);
                etamin = PHTelem{indexPatch}(i).vertex(2);
                etamax = PHTelem{indexPatch}(i).vertex(4);

                coordt = cell(ny,nx);
                dispmatx = zeros(ny,nx);
                dispmaty = zeros(ny,nx);
                stressvect = cell(ny,nx);
                stressVM = cell(ny,nx);
                wmatx = zeros(ny, nx);

                nument = size(PHTelem{indexPatch}(i).C,1);
                nodes = PHTelem{indexPatch}(i).nodes(1:nument);
                sctrx = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
                dsctrx = reshape([2*sctrx-1;2*sctrx],1,2*nument);
                cpts = controlPts{indexPatch}(nodes,1:2);
                wgts = controlPts{indexPatch}(nodes,3);
                for ii=1:nx
                    for jj=1:ny
                        dRdx = (PHTelem{indexPatch}(i).C)*squeeze(dBdu(ii,jj,:))*2/(ximax-ximin);
                        dRdy = (PHTelem{indexPatch}(i).C)*squeeze(dBdv(ii,jj,:))*2/(etamax-etamin);
                        RR = (PHTelem{indexPatch}(i).C)*squeeze(R(ii,jj,:));

                        RR = RR .* wgts;
                        dRdx = dRdx .* wgts;
                        dRdy = dRdy .* wgts;

                        w_sum = sum(RR);
                        dw_xi = sum(dRdx);
                        dw_eta = sum(dRdy);

                        dRdx = dRdx/w_sum - RR*dw_xi/w_sum^2;
                        dRdy = dRdy/w_sum - RR*dw_eta/w_sum^2;
                        RR = RR/w_sum;

                        % multiply by the jacobian of the transformation from reference
                        % space to the parameter space
                        coord = RR'*cpts;
                        dR  = [dRdx';dRdy'];
                        dxdxi = dR*cpts;


                        % Solve for first derivatives in global coordinates
                        dRdx = dxdxi\dR;

                        Bu = zeros(3,dim*nument);
                        for inode=1:nument
                            % Calculation of Bu
                            Bu(1,2*inode-1)=dRdx(1,inode);
                            Bu(2,2*inode)=dRdx(2,inode);
                            Bu(3,2*inode-1)=dRdx(2,inode);
                            Bu(3,2*inode)=dRdx(1,inode);
                        end

                        coordt{jj,ii} = coord;


                        % Calculate tdisp values
                        dispmatx(jj,ii) = dispmatx(jj,ii) + RR'*tdisp(2*sctrx-1);
                        dispmaty(jj,ii) = dispmaty(jj,ii) + RR'*tdisp(2*sctrx);

                        % Calculate stress values
                        stressvect{jj,ii} = material.C*Bu*tdisp(dsctrx);
                        stressVM{jj,ii} = sqrt(stressvect{jj,ii}(1)^2 - stressvect{jj,ii}(1)*stressvect{jj,ii}(2) + stressvect{jj,ii}(2)^2 +3*stressvect{jj,ii}(3)^2);

                        %calculate phase field values
                        wmatx(jj,ii) = wmatx(jj,ii) + RR'*tdisp(2*sizeBasis+sctrx);

                    end
                end
                physcoord((elementCounter-1)*4+1,:) = coordt{1,1};
                physcoord((elementCounter-1)*4+2,:) = coordt{1,2};
                physcoord((elementCounter-1)*4+3,:) = coordt{2,2};
                physcoord((elementCounter-1)*4+4,:) = coordt{2,1};

                dispcoord((elementCounter-1)*4+1,:) = [dispmatx(1,1) dispmaty(1,1)];
                dispcoord((elementCounter-1)*4+2,:) = [dispmatx(1,2) dispmaty(1,2)];
                dispcoord((elementCounter-1)*4+3,:) = [dispmatx(2,2) dispmaty(2,2)];
                dispcoord((elementCounter-1)*4+4,:) = [dispmatx(2,1) dispmaty(2,1)];


                sigmacoord((elementCounter-1)*4+1, 1) = stressVM{1,1};
                sigmacoord((elementCounter-1)*4+2, 1) = stressVM{1,2}';
                sigmacoord((elementCounter-1)*4+3, 1) = stressVM{2,2}';
                sigmacoord((elementCounter-1)*4+4, 1) = stressVM{2,1}';

                wcoord((elementCounter-1)*4+1) = wmatx(1,1);
                wcoord((elementCounter-1)*4+2) = wmatx(1,2);
                wcoord((elementCounter-1)*4+3) = wmatx(2,2);
                wcoord((elementCounter-1)*4+4) = wmatx(2,1);

            end
        end
    end

    % print onScreen or VTK
    if rem(globdat.ts.step,props.postProc.onScreen) == 0 || ...
           globdat.active == 0

        plot1 = subplot(2,2,1);
        cla(plot1)
        plotMesh(PHTelem,controlPts,geometry,0)
        axis equal

        plot2 = subplot(2,2,2);
        cla(plot2)
        plot(abs(lodi(1:step+1,2)),lodi(1:step+1,1),'-k','LineWidth',2.0)
        hold on;
        plot(abs(lodi(step+1,2)),lodi(step+1,1),'*r','MarkerSize',6.0)
        xlabel('Displacement [mm]')
        ylabel('Load [N]')
        xlim([0 PostProc.lodi.xmax])
        ylim([0 PostProc.lodi.ymax])
        set(gca,'FontSize',9)

        % plot3 = subplot(2,2,3);
        % cla(plot3)
        % trisurf(vIEN,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), log10(sigmacoord), 'EdgeColor','none','facecolor','interp')
        % view(0,90)
        % title('Displacements and \sigma_{VM}')
        % colorbar('vert')
        % colormap('parula')
        % hold on
        % axis tight
        % axis equal

        plot3 = subplot(2,2,3);
        cla(plot3)
        plot(errvec,'-k','LineWidth',2.0)
        hold on;
        xlabel('Iterations')
        ylabel('Error')
        set(gca, 'YScale', 'log')
        set(gca,'FontSize',9)

        plot4 = subplot(2,2,4);
        cla(plot4)
        factor = 0;
        trisurf(vIEN,physcoord(:,1)+dispcoord(:,1)*factor, physcoord(:,2)+dispcoord(:,2)*factor, zeros(size(physcoord,1),1), wcoord, 'facecolor','interp', 'EdgeColor','none')
        view(0,90)
        title('Phase Field')
        colorbar('vert')
        colormap('parula')
        hold on
        axis tight
        axis equal
        drawnow

        print(gcf,fullfile(dir_output,['Step',num2str(step),'.png']),'-dpng','-r300')

    end

    numVertices=4*size(vIEN,1);
    numVCtrlElmt  = 4;        
    writeVtk(numVertices,globdat.mesh.numElements,numVCtrlElmt,vIEN,physcoord,globdat.ts.step,wcoord,dispcoord,dir_output);           

end

end