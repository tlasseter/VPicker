

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

public class StsConsoleStatus implements StsStatusUI
{

    public StsConsoleStatus()
    {
    }

    public void setMaximum(float max) { }
    public void setMinimum(float min) { }
    public void setProgress(float n) { setText(String.valueOf(n)); }

    public void setMaximum(int max) { }
    public void setMinimum(int min) { }
    public void setProgress(int n) { setText(String.valueOf(n)); }
    public int getProgress() { return 0; }
    public void setTitle(String msg) { setText(msg);}
    public void setText(String msg) { System.out.println(msg);}
    public void setText(String msg, int msec) { setText(msg); sleep(msec);}
    public void sleep(int msec) { try { Thread.sleep(msec); } catch(Exception e) { }}
}
