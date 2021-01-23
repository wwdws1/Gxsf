## Tools说明  
- `DECAY.DEC`为总的Decay Card，自己写Decay Card时应先看此文件中是否有DIY产生子模型，若没有或为PHSP，则使用相空间产生MC样本，否则应使用已有的DIY产生子模型产生。   
- `pdt.table`中有各个例子对应的ID，在写MC Truth时需要参考。   
- `significance.cxx`为显著性计算工具，使用方式为`root ./significance.cxx`。\(Roofit得到的为Likelihood值，FDC得到的为s值\)  