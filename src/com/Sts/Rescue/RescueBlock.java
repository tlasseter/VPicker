package com.Sts.Rescue;
/*****************************************************************

    RescueBlock.java

    Rod Hanks   April 2000

******************************************************************/
public class RescueBlock
{
  public RescueBlock(String newBlockName,
									  RescueModel parentModel)
  {
    blockNdx = CreateBlock(newBlockName, parentModel.nativeNdx);
  }
  private native long CreateBlock(String newBlockName,
                                 long parentModel);
  protected long blockNdx = 0;
}
 
