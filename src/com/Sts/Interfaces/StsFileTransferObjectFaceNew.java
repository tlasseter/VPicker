package com.Sts.Interfaces;

/**
 * <p>Title: S2S development</p>
 * <p>Description: Integrated seismic to simulation software</p>
 * <p>Copyright: Copyright (c) 2005</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version c62e
 */
public interface StsFileTransferObjectFaceNew extends StsBasicFileTransferFaceNew
{
    public boolean hasReloadButton();
    public boolean hasOverrideButton();
    public boolean hasArchiveItButton();
	public void setReload(boolean reload);
	public boolean getReload();
    public void setArchiveIt(boolean archiveIt);
	public boolean getArchiveIt();
    public void setOverrideFilter(boolean override);
	public boolean getOverrideFilter();
}