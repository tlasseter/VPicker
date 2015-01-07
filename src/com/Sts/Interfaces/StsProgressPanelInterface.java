package com.Sts.Interfaces;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public interface StsProgressPanelInterface
{
   public void appendLine(String line);
   public void initializeProgressBar(int maxValue);
   public void setProgressBarMax(int maxValue);
   public void setProgressBarValue(int value);
}
