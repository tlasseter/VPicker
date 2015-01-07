package com.Sts.UI.Icons;

import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.image.*;

/**
 * Title:        Workflow development
 * Description:
 * Copyright:    Copyright (c) 2001
 * Company:      4D Systems LLC
 * @author T.J.Lasseter
 * @version 1.0
 */

public class StsIcon
{
    static StsIcon iconInstance = new StsIcon();
    static Class iconClass = iconInstance.getClass();

    public StsIcon()
    {
    }

    static public ImageIcon createIcon(String name)
	{
        return createIcon(iconClass, name);
    }
    
    static public Image createImage(String name)
	{
        return createImage(iconClass, name);
    }
    
    static public Image createImage(Class c, String filename)
    {
        try
        {
            java.net.URL url = c.getResource(filename);
            if(url == null) throw new StsException();
            Toolkit tk = Toolkit.getDefaultToolkit();
            return tk.createImage((ImageProducer)url.getContent());
        }
        catch (Exception e)
        {
            return null;
        }
    }
    static public ImageIcon createIcon(Class c, String filename)
    {
        try
        {
            Image img = createImage(c, filename);
            if(img == null) throw new StsException();
            return new ImageIcon(img);
        }
        catch (Exception e)
        {
//                  StsException.systemError("StsIcon.createIcon() failed.\n" +
//                          "Icon file not found: " + filename);
            return null;
        }
    }
}