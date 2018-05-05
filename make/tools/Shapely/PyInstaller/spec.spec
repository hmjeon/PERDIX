# -*- mode: python -*-

block_cipher = None


a = Analysis(['Shapely.py'],
             pathex=['C:\\Users\\Hyungmin\\Desktop\\Coding\\PERDIX-2L\\make\\tools\\Shapely\\pyinstaller'],
             binaries=[],
             datas=[],
             hiddenimports=['shapely'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
			 
a.binaries = a.binaries + [('geos_c.dll', 'geos_c.dll', 'BINARY')]
a.binaries = a.binaries + [('libgeos.lib', 'libgeos.lib', 'BINARY')]

pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='Shapely',
          debug=False,
          strip=False,
          upx=True,
          console=True )
