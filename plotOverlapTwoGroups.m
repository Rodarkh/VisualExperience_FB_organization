function [h, out, dataPoints]=plotOverlapTwoGroups(binMeans1, binMeans2, distantThresh,binAngle,col,yloc,g1,g2)

[H1, V1, M1]=getAngularAverages_dRF_Data(binMeans1, distantThresh,binAngle);
[H2, V2, M2]=getAngularAverages_dRF_Data(binMeans2, distantThresh,binAngle);


h=plot([1 2],[V1/(H1+V1) , V2/(H2+V2)],col)
ax=gca;
ylim(ax,[0.15 0.6])
xlim(ax,[0.9 2.1])
set(ax,'Xtick',1:2,'XtickLabel',{g1, g2})
ylabel(ax, 'Survival VertBin/(VertBin+HorzBin)')
text(ax, 1.5, yloc,['(',g1,') / (',g2,') =' num2str((V1/(H1+V1))/(V2/(H2+V2)))],'Color',col)

out = (V1/(H1+V1))/(V2/(H2+V2));
dataPoints = [V1/(H1+V1) , V2/(H2+V2)];
