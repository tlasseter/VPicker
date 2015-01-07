package com.Sts.Interfaces;

/**
 * Title:        Workflow development
 * Description:
 * Copyright:    Copyright (c) 2001
 * Company:      4D Systems LLC
 * @author T.J.Lasseter
 * @version 1.0
 */


public interface StsTimeEnabled
{
	public String getBornDate();
    public void setBornDate(String date);
    public String getDeathDate();
    public void setDeathDate(String date);
    public boolean isAlive(long projectTime);
}
