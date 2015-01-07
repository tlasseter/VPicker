package com.Sts.Actions.Wizards.AncillaryData;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;

import java.awt.*;


/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */
public class StsAncillaryDataDefine extends StsWizardStep
{
    StsAncillaryDataDefinePanel panel;
    StsHeaderPanel header;

    public StsAncillaryDataDefine(StsWizard wizard)
    {
        super(wizard);
        panel = new StsAncillaryDataDefinePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 400));
        header.setTitle("Ancillary Data File Selection");
        header.setSubtitle("Define File and Launch Instructions");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#AncillaryData");
        header.addInfoLine(wizardDialog, "(1) Select files from the table on the left.");
        header.addInfoLine(wizardDialog, "(2) Edit file definition and launch parameters");
        header.addInfoLine(wizardDialog, "(3) Once all file parameters are set press update button");
        header.addInfoLine(wizardDialog, "(4) Check the Assign to Well option if desired");
        header.addInfoLine(wizardDialog, "(5) Press the Next>> Button.");
    }

    public boolean start()
    {
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        return true;
    }

    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
  {
      if(panel.getSelectedObject() == null)
          return true;

      StsAncillaryData object = (StsAncillaryData)panel.getSelectedObject();
      object.defineLocation((StsGLPanel3d)glPanel, mouse);
      panel.updateBeans();
      return true;
  }

  public boolean assignToWell() { return panel.getAssignToWell(); }

}
