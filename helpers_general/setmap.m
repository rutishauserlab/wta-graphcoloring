function colorMapWeightMatrix=setmap()

load('colormap-weightMatrix-v2');
colorMapWeightMatrix=colormap2;

return;

%adjusted from redblue from colormap.org
colorMapWeightMatrix=[ 0.3711         0         0
    0.3771    0.0001    0.0001
    0.3831    0.0001    0.0001
    0.3890    0.0002    0.0002
    0.3950    0.0002    0.0002
    0.4010    0.0003    0.0003
    0.4070    0.0004    0.0004
    0.4130    0.0004    0.0004
    0.4189    0.0005    0.0005
    0.4249    0.0005    0.0005
    0.4309    0.0006    0.0006
    0.4369    0.0007    0.0007
    0.4429    0.0007    0.0007
    0.4489    0.0008    0.0008
    0.4548    0.0009    0.0009
    0.4608    0.0009    0.0009
    0.4668    0.0010    0.0010
    0.4728    0.0010    0.0010
    0.4788    0.0011    0.0011
    0.4847    0.0012    0.0012
    0.4907    0.0012    0.0012
    0.4967    0.0013    0.0013
    0.5027    0.0013    0.0013
    0.5087    0.0014    0.0014
    0.5146    0.0015    0.0015
    0.5206    0.0015    0.0015
    0.5266    0.0016    0.0016
    0.5326    0.0016    0.0016
    0.5386    0.0017    0.0017
    0.5446    0.0018    0.0018
    0.5505    0.0018    0.0018
    0.5565    0.0019    0.0019
    0.5625    0.0020    0.0020
    0.5685    0.0020    0.0020
    0.5745    0.0021    0.0021
    0.5804    0.0021    0.0021
    0.5864    0.0022    0.0022
    0.5924    0.0023    0.0023
    0.5984    0.0023    0.0023
    0.6044    0.0024    0.0024
    0.6104    0.0024    0.0024
    0.6163    0.0025    0.0025
    0.6223    0.0026    0.0026
    0.6283    0.0026    0.0026
    0.6343    0.0027    0.0027
    0.6403    0.0027    0.0027
    0.6462    0.0028    0.0028
    0.6522    0.0029    0.0029
    0.6582    0.0029    0.0029
    0.6642    0.0030    0.0030
    0.6702    0.0031    0.0031
    0.6761    0.0031    0.0031
    0.6821    0.0032    0.0032
    0.6881    0.0032    0.0032
    0.6941    0.0033    0.0033
    0.7001    0.0034    0.0034
    0.7061    0.0034    0.0034
    0.7120    0.0035    0.0035
    0.7180    0.0035    0.0035
    0.7240    0.0036    0.0036
    0.7300    0.0037    0.0037
    0.7360    0.0037    0.0037
    0.7419    0.0038    0.0038
    0.7479    0.0038    0.0038
    0.7539    0.0039    0.0039
    0.7559    0.0122    0.0122
    0.7580    0.0205    0.0205
    0.7600    0.0287    0.0287
    0.7620    0.0370    0.0370
    0.7641    0.0453    0.0453
    0.7661    0.0535    0.0535
    0.7682    0.0618    0.0618
    0.7702    0.0701    0.0701
    0.7722    0.0784    0.0784
    0.7743    0.0866    0.0866
    0.7763    0.0949    0.0949
    0.7783    0.1032    0.1032
    0.7804    0.1114    0.1114
    0.7824    0.1197    0.1197
    0.7844    0.1280    0.1280
    0.7865    0.1363    0.1363
    0.7885    0.1445    0.1445
    0.7905    0.1528    0.1528
    0.7926    0.1611    0.1611
    0.7946    0.1693    0.1693
    0.7966    0.1776    0.1776
    0.7987    0.1859    0.1859
    0.8007    0.1942    0.1942
    0.8028    0.2024    0.2024
    0.8048    0.2107    0.2107
    0.8068    0.2190    0.2190
    0.8089    0.2273    0.2273
    0.8109    0.2355    0.2355
    0.8129    0.2438    0.2438
    0.8150    0.2521    0.2521
    0.8170    0.2603    0.2603
    0.8190    0.2686    0.2686
    0.8211    0.2769    0.2769
    0.8231    0.2852    0.2852
    0.8251    0.2934    0.2934
    0.8272    0.3017    0.3017
    0.8292    0.3100    0.3100
    0.8312    0.3182    0.3182
    0.8333    0.3265    0.3265
    0.8353    0.3348    0.3348
    0.8373    0.3431    0.3431
    0.8394    0.3513    0.3513
    0.8414    0.3596    0.3596
    0.8435    0.3679    0.3679
    0.8455    0.3761    0.3761
    0.8475    0.3844    0.3844
    0.8496    0.3927    0.3927
    0.8516    0.4010    0.4010
    0.8536    0.4092    0.4092
    0.8557    0.4175    0.4175
    0.8577    0.4258    0.4258
    0.8597    0.4341    0.4341
    0.8618    0.4423    0.4423
    0.8638    0.4506    0.4506
    0.8658    0.4589    0.4589
    0.8679    0.4671    0.4671
    0.8699    0.4754    0.4754
    0.8719    0.4837    0.4837
    0.8740    0.4920    0.4920
    0.8760    0.5002    0.5002
    0.8781    0.5085    0.5085
    0.8801    0.5168    0.5168
    0.8821    0.5250    0.5250
    0.8842    0.5333    0.5333
    0.8862    0.5416    0.5416
    0.8882    0.5499    0.5499
    0.8903    0.5581    0.5581
    0.8923    0.5664    0.5664
    0.8943    0.5747    0.5747
    0.8964    0.5830    0.5830
    0.8984    0.5912    0.5912
    0.9004    0.5995    0.5995
    0.9025    0.6078    0.6078
    0.9045    0.6160    0.6160
    0.9065    0.6243    0.6243
    0.9086    0.6326    0.6326
    0.9106    0.6409    0.6409
    0.9127    0.6491    0.6491
    0.9147    0.6574    0.6574
    0.9167    0.6657    0.6657
    0.9188    0.6739    0.6739
    0.9208    0.6822    0.6822
    0.9228    0.6905    0.6905
    0.9249    0.6988    0.6988
    0.9269    0.7070    0.7070
    0.9289    0.7153    0.7153
    0.9310    0.7236    0.7236
    0.9330    0.7318    0.7318
    0.9350    0.7401    0.7401
    0.9371    0.7484    0.7484
    0.9391    0.7567    0.7567
    0.9411    0.7649    0.7649
    0.9432    0.7732    0.7732
    0.9452    0.7815    0.7815
    0.9472    0.7898    0.7898
    0.9493    0.7980    0.7980
    0.9513    0.8063    0.8063
    0.9534    0.8146    0.8146
    0.9554    0.8228    0.8228
    0.9574    0.8311    0.8311
    0.9595    0.8394    0.8394
    0.9615    0.8477    0.8477
    0.9635    0.8559    0.8559
    0.9656    0.8642    0.8642
    0.9676    0.8725    0.8725
    0.9696    0.8807    0.8807
    0.9717    0.8890    0.8890
    0.9737    0.8973    0.8973
    0.9757    0.9056    0.9056
    0.9778    0.9138    0.9138
    0.9798    0.9221    0.9221
    0.9818    0.9304    0.9304
    0.9839    0.9386    0.9386
    0.9859    0.9469    0.9469
    0.9880    0.9552    0.9552
    0.9900    0.9635    0.9635
    0.9920    0.9717    0.9717
    0.9941    0.9800    0.9800
    0.9961    0.9883    0.9883
    0.9945    0.9883    0.9898
    0.9930    0.9883    0.9914
    0.9914    0.9883    0.9930
    0.9898    0.9883    0.9945
    0.9883    0.9883    0.9961
    0.9391    0.9391    0.9844
    0.8898    0.8898    0.9727
    0.8406    0.8406    0.9609
    0.7914    0.7914    0.9492
    0.7422    0.7422    0.9375
    0.6930    0.6930    0.9258
    0.6438    0.6438    0.9141
    0.5945    0.5945    0.9023
    0.5453    0.5453    0.8906
    0.4961    0.4961    0.8789
    0.4469    0.4469    0.8672
    0.3977    0.3977    0.8555
    0.3484    0.3484    0.8438
    0.2992    0.2992    0.8320
    0.2500    0.2500    0.8203
    0.2008    0.2008    0.8086
    0.1516    0.1516    0.7969
    0.1023    0.1023    0.7852
    0.0531    0.0531    0.7734
    0.0039    0.0039    0.7617
    0.0038    0.0038    0.7534
    0.0037    0.0037    0.7451
    0.0037    0.0037    0.7368
    0.0036    0.0036    0.7285
    0.0035    0.0035    0.7202
    0.0034    0.0034    0.7119
    0.0033    0.0033    0.7035
    0.0032    0.0032    0.6952
    0.0032    0.0032    0.6869
    0.0031    0.0031    0.6786
    0.0030    0.0030    0.6703
    0.0029    0.0029    0.6620
    0.0028    0.0028    0.6537
    0.0027    0.0027    0.6454
    0.0027    0.0027    0.6371
    0.0026    0.0026    0.6287
    0.0025    0.0025    0.6204
    0.0024    0.0024    0.6121
    0.0023    0.0023    0.6038
    0.0022    0.0022    0.5955
    0.0022    0.0022    0.5872
    0.0021    0.0021    0.5789
    0.0020    0.0020    0.5706
    0.0019    0.0019    0.5623
    0.0018    0.0018    0.5539
    0.0017    0.0017    0.5456
    0.0017    0.0017    0.5373
    0.0016    0.0016    0.5290
    0.0015    0.0015    0.5207
    0.0014    0.0014    0.5124
    0.0013    0.0013    0.5041
    0.0012    0.0012    0.4958
    0.0012    0.0012    0.4875
    0.0011    0.0011    0.4791
    0.0010    0.0010    0.4708
    0.0009    0.0009    0.4625
    0.0008    0.0008    0.4542
    0.0007    0.0007    0.4459
    0.0007    0.0007    0.4376
    0.0006    0.0006    0.4293
    0.0005    0.0005    0.4210
    0.0004    0.0004    0.4126
    0.0003    0.0003    0.4043
    0.0002    0.0002    0.3960
    0.0002    0.0002    0.3877
    0.0001    0.0001    0.3794
         0         0    0.3711];
     