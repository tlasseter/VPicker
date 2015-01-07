package com.Sts.Workflow;

import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Icons.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsSplashPanel extends StsJPanel
{
    public StsSplashPanel()
    {
        setBackground(Color.BLACK);
        JLabel iconLabel = new JLabel(StsIcon.createIcon(Main.vendorName + "SplashScreen.jpg"));
        add(iconLabel);
    }
}