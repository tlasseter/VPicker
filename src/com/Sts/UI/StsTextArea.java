package com.Sts.UI;

import javax.swing.*;
import java.awt.*;

/**
 * Title:        Workflow development
 * Description:
 * Copyright:    Copyright (c) 2001
 * Company:      4D Systems LLC
 * @author T.J.Lasseter
 * @version 1.0
 */

public class StsTextArea extends JTextArea
{
    public StsTextArea()
	{
        setBackground(Color.lightGray);
        setEditable(false);
        setText("");
        setLineWrap(false);
        setFont(new Font("Monospaced", 1, 10));
        setBorder(BorderFactory.createEtchedBorder());
        setDoubleBuffered(true);
        setCaretColor(Color.red);
        setCaretPosition(0);
    }
}