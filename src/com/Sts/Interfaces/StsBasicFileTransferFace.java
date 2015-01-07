package com.Sts.Interfaces;

import com.Sts.IO.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jun 2, 2008
 * Time: 9:43:54 AM
 * To change this template use File | Settings | File Templates.
 */
public interface StsBasicFileTransferFace
{
    void addFiles(StsFile[] files);
    void removeFiles(StsFile[] files);
    void removeAllFiles();
    void fileSelected(StsFile selectedFile);
    void availableFileSelected(StsFile availableFile);
    boolean hasDirectorySelection();
}
