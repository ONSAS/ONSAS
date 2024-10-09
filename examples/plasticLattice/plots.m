load 'data.mat'

deltas  = -matUs((loadednode-1)*6+3,1:10:end)' ;
lambdas = loadFactorsMat(1:10:end,3) ;

figure
plot( deltas , lambdas, 'b-x' )
grid on, hold on
# plot( deltas , loadFactorsMat(:,2), 'r-o' )
# plot( deltasB , loadFactorsMatB(:,2), 'g-o' )
# plot( deltasB , valsPB, 'k-*' )
xlabel('displacement')
ylabel('load factor \lambda')
axis tight
print('fig.png')
# legend('analytic-hard','numeric-hard', 'analytic-soft','numeric-soft')
