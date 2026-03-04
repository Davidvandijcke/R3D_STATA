/*
 * diagnose_plugin.do - Diagnose plugin loading issues
 */

version 14.0
clear all

di as text _n "{hline 60}"
di as text "R3D Plugin Diagnostics"
di as text "{hline 60}"

* System information
di as text _n "System Information:"
di as text "  OS: " as result c(os)
di as text "  Machine type: " as result c(machine_type)
di as text "  Bit: " as result c(bit)
di as text "  Stata version: " as result c(stata_version)

* Check plugin file
di as text _n "Checking plugin file..."
local pluginpath "ado/r3d_plugin.plugin"
capture confirm file "`pluginpath'"
if _rc == 0 {
    di as result "✓ Plugin file found: `pluginpath'"
}
else {
    di as error "✗ Plugin file not found at: `pluginpath'"
    exit 1
}

* Try to load the plugin directly
di as text _n "Attempting to load plugin..."

* Method 1: Direct load
capture program r3d_plugin, plugin
if _rc == 0 {
    di as result "✓ Plugin loaded successfully with direct method"
    program drop r3d_plugin
}
else {
    di as text "⚠ Direct load failed (rc=" _rc ")"
    
    * Method 2: Load with full path
    local fullpath = c(pwd) + "/ado/r3d_plugin.plugin"
    di as text "  Trying with full path: " as result "`fullpath'"
    capture program r3d_plugin, plugin using("`fullpath'")
    if _rc == 0 {
        di as result "✓ Plugin loaded with full path"
        program drop r3d_plugin
    }
    else {
        di as error "✗ Plugin load failed with full path (rc=" _rc ")"
    }
}

* Load the plugin loader
di as text _n "Testing plugin loader..."
quietly capture do ado/r3d_plugin_load.ado
if _rc == 0 {
    di as result "✓ Plugin loader executed"
    
    if "$R3D_USE_PLUGIN" == "1" {
        di as result "✓ R3D_USE_PLUGIN = 1 (Plugin available)"
    }
    else {
        di as text "⚠ R3D_USE_PLUGIN = 0 (Plugin not available)"
    }
}
else {
    di as error "✗ Failed to load plugin loader"
}

di as text _n "{hline 60}"