                (
                    traction[patchI]

                  - (n &
                    (
                        mu.value()*gradUpatches[patchI].T()
                      - (mu + lambda).value()*gradUpatches[patchI]
                    )
                    )

                  - n*
                    (
                        lambda.value()*tr(gradUpatches[patchI])

                    )

                )/(2.0*mu + lambda).value();
