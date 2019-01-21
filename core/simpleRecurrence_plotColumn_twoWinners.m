function simpleRecurrence_plotColumn_twoWinners(till, x1, y1, inputAmp, subplotNrs,pstr)

ax1=subplot(3,3,subplotNrs(1));
x=1:till;
plot( x, x1(1,:), 'r', x, x1(2,:), 'b',x, y1(1,:), 'c' );
title([pstr]);
legend('x1','x2','inh');

ax2=subplot(3,3,subplotNrs(2));
x=1:till;
plot( x(2:end), diff(x1(1,:)), 'r', x(2:end), diff(x1(2,:)), 'b',  x(2:end), diff(y1(1,:)), 'c' );
title(['diff ' pstr ]);
legend('x1','x2','inh');

ax3=subplot(3,3,subplotNrs(3));
x=1:till;
d1=x1(1,:)/inputAmp(1);
d2=x1(2,:)/inputAmp(2);

plot( x, d1, 'r', x, d2, 'b');
%ylim([-0.1 0.1]);
ylabel('gain');
legend('x1','x2');



linkaxes([ax1 ax2 ax3],'x');