for i=1:nmorphs
    spikepeak1(i) = floor(bval(i).info(1).actualpeaks(1) * 24000) ;
end
figure ;
colours = 'rgbcmrgbcmrgbcm' ;
for i=1:nmorphs
    plot(bval(i).info(1).final(spikepeak1(i)-50:spikepeak1(i)+44),colours(i)) ; hold on ;
end
spikepeak2 = floor(bval(1).info(2).actualpeaks(1) * 24000) ;
plot(bval(nmorphs).info(2).final(spikepeak2-50:spikepeak2+44),'k') ;

