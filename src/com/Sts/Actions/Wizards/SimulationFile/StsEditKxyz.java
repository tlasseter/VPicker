package com.Sts.Actions.Wizards.SimulationFile;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.RFUI.*;

public class StsEditKxyz extends StsWizardStep
{
    StsKxyzPanel panel;
    StsHeaderPanel header;

    public StsEditKxyz(StsWizard wizard)
    {
        super(wizard, new StsKxyzPanel(), null, new StsHeaderPanel());
        panel = (StsKxyzPanel) getContainer();
        panel.setKx(1.0);
        panel.setKy(1.0);
        panel.setKz(0.1);
//        panel.setTitle("Edit the variogram parameters:");
//		panel.setPreferredSize(new Dimension(400,300));
        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Simulation File Export");
        header.setSubtitle("Edit Kxyz");
    }

    public boolean start()
    {
        return true;
    }

    public boolean end()
    {
        return true;
    }

    public double getKx() { return panel.getKx(); }
    public double getKy() { return panel.getKy(); }
    public double getKz() { return panel.getKz(); }

    public void setKx(double kx) { panel.setKx(kx); }
    public void setKy(double ky) { panel.setKx(ky); }
    public void setKz(double kz) { panel.setKx(kz); }
}
