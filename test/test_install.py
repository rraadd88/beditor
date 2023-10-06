def test_install():
    try:
        from beditor import pipeline
        # cfg=get_deps(cfg)
        # cfg=get_genomes(cfg)
        print(f">>> SUCCESS")
    except:
        print(f">>> TEST NOT SUCCESSFUL. Something's wrong with installation.")

test_install()