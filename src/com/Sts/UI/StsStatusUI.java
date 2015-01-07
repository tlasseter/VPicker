
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;



public interface StsStatusUI
{
    public void setMaximum(float max);
    public void setMinimum(float min);
    public void setProgress(float n);

    public void setMaximum(int max);
    public void setMinimum(int min);
    public void setProgress(int n);
    public int getProgress();
    public void setTitle(String msg);
    public void setText(String msg);
    public void setText(String msg, int msec);
    public void sleep(int msec);
}
