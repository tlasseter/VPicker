package com.Sts.Interfaces;

/**
 * Title:
 * Description:
 * Copyright:    Copyright (c) 2001
 * Company:
 * @author
 * @version 1.0
 */

import javax.swing.*;

public interface StsListModel extends ListModel
{
	String getLongestString();
	void addLine(String line);
}