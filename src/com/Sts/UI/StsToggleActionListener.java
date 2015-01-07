package com.Sts.UI;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

import com.Sts.Reflect.*;
import com.Sts.Utilities.*;

import java.awt.event.*;
import java.lang.reflect.*;

public class StsToggleActionListener implements ActionListener
{
    protected StsToggleButton toggleButton;
  	protected StsMethod selectMethod, deselectMethod;

    static final boolean debug = false;

    public StsToggleActionListener(StsToggleButton toggleButton, StsMethod selectMethod, StsMethod deselectMethod)
    {
        this.toggleButton = toggleButton;
    	this.selectMethod = selectMethod;
        this.deselectMethod = deselectMethod;
    }

    public void actionPerformed(ActionEvent event)
    {
        StsMethod useMethod;

        if(toggleButton.getModel().isSelected())
            useMethod = selectMethod;
        else
            useMethod = deselectMethod;

        try
        {
            if(useMethod == null) return;
            Method method = useMethod.getMethod();
            if(method == null) return;
            useMethod.getMethod().invoke(useMethod.getInstance(), useMethod.getMethodArgs());
            if(debug) System.out.println("StsToggleActionListener.actionPerformed() toggle.isSelected is: " + toggleButton.isSelected());
        }
        catch (Exception e)
        {
            StsException.outputException( "method=" + useMethod.getMethod().getName(), e, StsException.WARNING);
        }
    }
}
