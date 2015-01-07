package com.Sts.Types;

import com.Sts.MVC.ViewWell.*;
import com.Sts.MVC.*;
import com.Sts.DBTypes.*;
import com.Sts.Actions.Wizards.*;

import javax.swing.*;
import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Dec 17, 2009
 * Time: 8:38:27 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsWellWizardViewModel extends StsWellViewModel
{
    /** wizard this wellModel and view belongs to; used when wellWindow is resized */
    transient StsWizard wizard = null;

    public Frame getParentFrame() { return null; }

    protected Component getCenterComponent() { return wizard.dialog; }

    public StsWellWizardViewModel()
    {
    }

    public StsWellWizardViewModel(StsWell well, StsWizard wizard)
    {
        super(well, wizard.model);
        this.wizard = wizard;
        initializeWellWindowPanel();
    }

    public void rebuild()
    {
        wizard.rebuild();
    }
}